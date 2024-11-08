#!/usr/bin/env python3
"""
Copyright (c) 2024 Zhejiang Lab

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
"""

import multiprocessing
import threading
import queue
import logging
import time
import os
import signal
import random
from queue import Empty
from ctypes import *
from functools import wraps


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s')

# Load dynamic link library
lib_path = os.getenv('LD_LIBRARY_PATH', '/home/soft/presto/lib')
lib_file = os.path.join(lib_path, 'libaccelsearch_lib.so')

if not os.path.exists(lib_file):
    raise FileNotFoundError(f"Shared library not found: {lib_file}")

lib = CDLL(lib_file)


# Define fcomplex type
class FComplex(Structure):
    _fields_ = [("r", c_float),
                ("i", c_float)]

# Define Kernel structure
class Kernel(Structure):
    _fields_ = [
        ("z", c_int),
        ("w", c_int),
        ("fftlen", c_int),
        ("numgoodbins", c_int),
        ("numbetween", c_int),
        ("kern_half_width", c_int),
        ("data", POINTER(FComplex))
    ]

# Define Subharminfo structure
class Subharminfo(Structure):
    _fields_ = [
        ("numharm", c_int),
        ("harmnum", c_int),
        ("zmax", c_int),
        ("wmax", c_int),
        ("numkern_zdim", c_int),
        ("numkern_wdim", c_int),
        ("numkern", c_int),
        ("kern", POINTER(POINTER(Kernel))),
        ("rinds", POINTER(c_ushort)),
        ("zinds", POINTER(c_ushort))
    ]

# Define accelobs structure
class AccelObs(Structure):
    _fields_ = [
        ("N", c_longlong),
        ("numbins", c_longlong),
        ("lobin", c_longlong),
        ("highestbin", c_longlong),
        ("maxkernlen", c_int),
        ("corr_uselen", c_int),
        ("fftlen", c_int),
        ("numharmstages", c_int),
        ("numz", c_int),
        ("numw", c_int),
        ("numbetween", c_int),
        ("numzap", c_int),
        ("dat_input", c_int),
        ("mmap_file", c_int),
        ("inmem", c_int),
        ("norm_type", c_int),
        ("dt", c_double),
        ("T", c_double),
        ("rlo", c_double),
        ("rhi", c_double),
        ("dr", c_double),
        ("zlo", c_double),
        ("zhi", c_double),
        ("dz", c_double),
        ("wlo", c_double),
        ("whi", c_double),
        ("dw", c_double),
        ("baryv", c_double),
        ("nph", c_float),
        ("sigma", c_float),
        ("powcut", POINTER(c_float)),
        ("ffdotplane", POINTER(c_float)),
        ("lobins", POINTER(c_double)),
        ("hibins", POINTER(c_double)),
        ("numindep", POINTER(c_longlong)),
        ("fftfile", c_void_p),
        ("workfile", c_void_p),
        ("fft", POINTER(FComplex)),
        ("rootfilenm", c_char_p),
        ("candnm", c_char_p),
        ("accelnm", c_char_p),
        ("workfilenm", c_char_p),
        ("use_harmonic_polishing", c_int)
    ]

# Define Infodata structure
class Infodata(Structure):
    _fields_ = [
        ("ra_s", c_double),
        ("dec_s", c_double),
        ("N", c_double),
        ("dt", c_double),
        ("fov", c_double),
        ("mjd_f", c_double),
        ("dm", c_double),
        ("freq", c_double),
        ("freqband", c_double),
        ("chan_wid", c_double),
        ("wavelen", c_double),
        ("waveband", c_double),
        ("energy", c_double),
        ("energyband", c_double),
        ("onoff", c_double * 80),  
        ("num_chan", c_int),
        ("mjd_i", c_int),
        ("ra_h", c_int),
        ("ra_m", c_int),
        ("dec_d", c_int),
        ("dec_m", c_int),
        ("bary", c_int),
        ("numonoff", c_int),
        ("notes", c_char * 500),
        ("name", c_char * 200),
        ("object", c_char * 100),
        ("instrument", c_char * 100),
        ("observer", c_char * 100),
        ("analyzer", c_char * 100),
        ("telescope", c_char * 40),
        ("band", c_char * 40),
        ("filt", c_char * 7)
    ]

# Complete Cmdline structure definition
class Cmdline(Structure):
    _fields_ = [
        ("batchsizeP", c_char),
        ("batchsize", c_int),
        ("batchsizeC", c_int),
        ("ncpusP", c_char),
        ("ncpus", c_int),
        ("ncpusC", c_int),
        ("lobinP", c_char),
        ("lobin", c_int),
        ("lobinC", c_int),
        ("numharmP", c_char),
        ("numharm", c_int),
        ("numharmC", c_int),
        ("zmaxP", c_char),
        ("zmax", c_int),
        ("zmaxC", c_int),
        ("wmaxP", c_char),
        ("wmax", c_int),
        ("wmaxC", c_int),
        ("sigmaP", c_char),
        ("sigma", c_float),
        ("sigmaC", c_int),
        ("rloP", c_char),
        ("rlo", c_double),
        ("rloC", c_int),
        ("rhiP", c_char),
        ("rhi", c_double),
        ("rhiC", c_int),
        ("floP", c_char),
        ("flo", c_double),
        ("floC", c_int),
        ("fhiP", c_char),
        ("fhi", c_double),
        ("fhiC", c_int),
        ("inmemP", c_char),
        ("photonP", c_char),
        ("medianP", c_char),
        ("locpowP", c_char),
        ("zaplistP", c_char),
        ("zaplist", c_char_p),
        ("zaplistC", c_int),
        ("baryvP", c_char),
        ("baryv", c_double),
        ("baryvC", c_int),
        ("otheroptP", c_char),
        ("noharmpolishP", c_char),
        ("noharmremoveP", c_char),
        ("argc", c_int),
        ("argv", POINTER(c_char_p)),
        ("full_cmd_line", c_char_p),
        ("cudaP", c_char),
        ("cuda", c_char)
    ]


# Define GSList structure
class GSList(Structure):
    pass

# GSList structure fields need to reference itself, so set them separately after class definition
GSList._fields_ = [
    ("data", c_void_p),   # gpointer represented by c_void_p
    ("next", POINTER(GSList))  # Pointer to the next GSList node
]

def random_sleep():
    # Generate a random integer from 0 to 2 (inclusive)
    if random.randint(0, 2) == 0:  # 1/3 chance to generate 0
        print("Sleeping for 20 seconds...")
        time.sleep(20)
    else:
        print("Not sleeping.")

def redirect_output_if_debug(debug_mode):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if not debug_mode:
                stdout = os.dup(1)  # Save standard output
                stderr = os.dup(2)  # Save standard error
                silent = os.open(os.devnull, os.O_WRONLY)  # Open /dev/null device

                os.dup2(silent, 1)  # Redirect standard output to /dev/null
                os.dup2(silent, 2)  # Redirect standard error to /dev/null

            result = func(*args, **kwargs)  # Call the original function

            if not debug_mode:
                os.dup2(stdout, 1)  # Restore standard output
                os.dup2(stderr, 2)  # Restore standard error
                os.close(silent)  # Close /dev/null device
                os.close(stdout)  # Close saved standard output
                os.close(stderr)  # Close saved standard error
            return result
        return wrapper
    return decorator


class Worker:
    def __init__(self, task_queue, task_status_dict, task_map, result_dict, gpu_id=0):
        self.task_queue = task_queue
        self.task_status_dict = task_status_dict
        self.task_map = task_map
        self.result_dict = result_dict
        self.stop_event = multiprocessing.Event()
        self.gpu_id = gpu_id
        self.debug_mode = os.getenv('DEBUG_MODE', 'false').lower() == 'true'     
        # Define function pointers and apply decorators
        self.accelsearch_CPU1 = redirect_output_if_debug(self.debug_mode)(CFUNCTYPE(None, c_int, POINTER(c_char_p), POINTER(POINTER(POINTER(Subharminfo))), POINTER(AccelObs), POINTER(Infodata), POINTER(POINTER(Cmdline)), use_errno=False)(('accelsearch_CPU1', lib)))
        self.accelsearch_GPU = redirect_output_if_debug(self.debug_mode)(CFUNCTYPE(c_int, AccelObs, POINTER(POINTER(Subharminfo)), POINTER(POINTER(GSList)), POINTER(Cmdline), use_errno=False)(('accelsearch_GPU', lib)))
        self.accelsearch_CPU2 = redirect_output_if_debug(self.debug_mode)(CFUNCTYPE(None, POINTER(POINTER(GSList)), POINTER(AccelObs), POINTER(Infodata), POINTER(Cmdline), use_errno=False)(('accelsearch_CPU2', lib)))                                                                
        self.create_time = time.time()
        self.timeout_event = multiprocessing.Event()  # Each worker has its own timeout event
        self.last_active_times = multiprocessing.Manager().dict()
        self.global_lock = multiprocessing.Lock()  # Global lock
        self.process = multiprocessing.Process(target=self.process_tasks)
        self.process.start()

    def update_timestamp(self, thread_id):
        with self.global_lock:
            self.last_active_times[thread_id] = time.time()

    def cancel_timestamp(self, thread_id):
        with self.global_lock:
            self.last_active_times[thread_id] = None

    def process_tasks(self):                                                                                                                 
        t1_t2_queue = queue.Queue(maxsize=1)
        t2_t3_queue = queue.Queue(maxsize=1)

        t1 = threading.Thread(target=self.task_t1, args=(self.task_queue, t1_t2_queue, 't1'))
        t2 = threading.Thread(target=self.task_t2, args=(t1_t2_queue, t2_t3_queue, self.result_dict, 't2'))
        t3 = threading.Thread(target=self.task_t3, args=(t2_t3_queue, self.result_dict, 't3'))

        t1.start()
        t2.start()
        t3.start()

        t1.join()
        t2.join()
        t3.join()

    def task_t1(self, input_queue, output_queue, thread_id):
        while not self.stop_event.is_set():
            try:
                data = input_queue.get(block=False, timeout=0.1)
                logging.debug(f"t1 get data: {data}")
                task_id, args = data
                with self.global_lock:
                    self.task_status_dict[task_id] = self.process.pid
                    self.task_map[task_id] = args

                start_time = time.time()  # Start timing
                logging.debug(f'Task T1 starts on {args[-1]}.')
                argc = len(args)
                argv = (c_char_p * argc)(*map(lambda s: s.encode('utf-8'), args))

                obs = AccelObs()
                idata = Infodata()
                cmd = Cmdline()
                subharminfs = POINTER(POINTER(Subharminfo))()
                pcmd = pointer(cmd)

                self.update_timestamp(thread_id)

                self.accelsearch_CPU1(argc, argv, byref(subharminfs), byref(obs), byref(idata), byref(pcmd))
                self.cancel_timestamp(thread_id)
                elapsed_time = time.time() - start_time  # Calculate elapsed time
                logging.debug(f'Task T1 executed in {elapsed_time:.2f} seconds.')  # Log elapsed time
                output_queue.put((task_id, subharminfs, obs, idata, pcmd, args))
            except Empty:
                 time.sleep(0.1)
            except Exception as e:
                logging.error(f"Error in task T1: {str(e)}")
        logging.info(f"Finish t1 in thread {self.process.name}")

    def task_t2(self, input_queue, output_queue, result_dict, thread_id):
        os.environ['CUDA_VISIBLE_DEVICES'] = str(self.gpu_id)
        while not self.stop_event.is_set():
            try:
                data = input_queue.get(block=False, timeout=0.1)
                logging.debug(f"t2 get data: {data}")
                start_time = time.time()  # Start timing
                task_id, subharminfs, obs, idata, pcmd, args = data
                logging.debug(f'Task T2 starts on {args[-1]}.')
                cands = POINTER(GSList)()
                self.update_timestamp(thread_id)
                error_code = self.accelsearch_GPU(obs, subharminfs, byref(cands), pcmd)
                # random_sleep()
                self.cancel_timestamp(thread_id)
                if error_code != 0:
                    logging.error(f"Task T2 completed with error_code: {error_code}, task_id: {task_id}, args: {args}")

                    with self.global_lock:
                        result_dict[task_id] = error_code
                        # Remove the task from task status dictionary once done
                        del self.task_status_dict[task_id]
                        del self.task_map[task_id]
                    continue

                elapsed_time = time.time() - start_time  # Calculate elapsed time
                logging.debug(f'Task T2 executed in {elapsed_time:.2f} seconds.')  # Log elapsed time
                output_queue.put((task_id, cands, obs, idata, pcmd, args))
            except Empty:
                 time.sleep(0.1)
            except Exception as e:
                logging.error(f"Error in task T2: {str(e)}")
        logging.info(f"Finish t2 in thread {self.process.name}")


    def task_t3(self, input_queue, result_dict, thread_id):
        while not self.stop_event.is_set():
            try:
                data = input_queue.get(block=False, timeout=0.1)
                logging.debug(f"t3 get data: {data}")
                start_time = time.time()  # Start timing
                task_id, cands, obs, idata, pcmd, args = data
                logging.debug(f'Task T3 starts on {args[-1]}.')
                self.update_timestamp(thread_id)
                self.accelsearch_CPU2(byref(cands), byref(obs), byref(idata), pcmd)
                self.cancel_timestamp(thread_id)

                with self.global_lock:
                    result_dict[task_id] = 0
                    # Remove the task from task status dictionary once done
                    del self.task_status_dict[task_id]
                    del self.task_map[task_id]
                logging.debug(f"result_dict in t3, task_id: {task_id}, result: {result_dict[task_id]}")
                elapsed_time = time.time() - start_time  # Calculate elapsed time
                logging.debug(f'Task T3 executed in {elapsed_time:.2f} seconds.')  # Log elapsed time
            except Empty:
                 time.sleep(0.1)
            except Exception as e:
                logging.error(f"Error in task T3: {str(e)}")
        logging.info(f"Finish t3 in thread {self.process.name}")

class Pool:
    def __init__(self, num_workers_per_gpu, timeout_seconds=30, recreate_pool_timeout_seconds=60*60, gpu_count=1):
        """
        Initialize the thread pool.

        :param num_workers_per_gpu: Number of worker threads per GPU in the pool.
        :param timeout_seconds: Task timeout in seconds.
        :param recreate_pool_timeout_seconds: Time interval to recreate the thread pool (seconds). Purpose: To avoid memory leaks.
        """
        manager = multiprocessing.Manager()
        self.task_queue = multiprocessing.Queue()
        self.task_status_dict = manager.dict()
        self.result_dict = manager.dict()
        self.task_map = manager.dict()
        self.stop_event = threading.Event()  # Event to signal stop
        self.timeout_seconds = timeout_seconds
        self.recreate_pool_timeout_seconds = recreate_pool_timeout_seconds
        self.workers = []
        # Create workers and assign them to the corresponding GPU
        for gpu_id in range(int(gpu_count)):
            for _ in range(int(num_workers_per_gpu)):
                self.workers.append(self.create_worker(gpu_id=gpu_id))
        self.monitor_thread = threading.Thread(target=self.monitor_workers)
        self.monitor_thread.daemon = True 
        self.monitor_thread.start()
        self.task_id_counter = 0
        self.lock = threading.Lock()

    def need_recreate_worker(self, worker):
        """
        Determine if a worker needs to be recreated, not thread-safe
        """
        current_time = time.time()
        # self.recreate_pool_timeout_seconds randomly fluctuates
        if current_time - worker.create_time > self.recreate_pool_timeout_seconds * random.uniform(0.8, 1.5):
            logging.error(f"Pool Timeout detected for Worker#{worker.process.name}.")
            return True
        for thread_id, last_active_time in list(worker.last_active_times.items()):
            if last_active_time and current_time - last_active_time > self.timeout_seconds:
                logging.error(f"Process Timeout detected for thread {thread_id} in Worker#{worker.process.name}.")
                return True
        return False

    def monitor_workers(self):
        logging.info("monitor workers...")
        while not self.stop_event.is_set():
            time.sleep(0.5)
            try:
                for i, worker in enumerate(self.workers):
                    if self.need_recreate_worker(worker):
                        with worker.global_lock:
                            with self.lock:
                                 # Since we are operating on workers[i], we need to lock

                                if worker.stop_event.is_set():
                                    break
                                wpid = worker.process.pid
                                worker.stop_event.set()
                                # sleep 10s and logging
                                logging.info(f"Waiting Worker#{worker.process.name}...")
                                timeout = 15  # Set timeout to 5 seconds
                                worker.process.join(timeout)
                                if worker.process.exitcode is None:
                                    logging.info(f"Terminating Worker#{worker.process.name}...")
                                    os.kill(worker.process.pid, signal.SIGKILL)  # Force kill the process
                                    logging.info(f"Joining Worker#{worker.process.name}...")
                                    worker.process.join()  # Wait for the process to terminate
                                    logging.info(f"Finish Worker#{worker.process.name}...")
                                for task_id, pid in list(self.task_status_dict.items()):
                                    if pid != wpid:
                                        continue
                                    task = self.task_map[task_id]
                                    logging.info(f"reschedule {(task_id, task)}")
                                    self.task_queue.put((task_id, task))
                                    del self.task_status_dict[task_id]
                                    del self.task_map[task_id]
                                self.workers[i] = self.create_worker(worker.gpu_id)  # Replace with a new worker
                            logging.info(f"Worker#{worker.process.name} is replaced with Worker#{self.workers[i].process.name}")
                            break
            except Exception as e:
                logging.error(f"An exception occurred in monitor_workers: {e}")
                continue  # Continue the while loop despite the exception

    def submit(self, cmd_args):
        self.task_id_counter += 1
        self.task_queue.put((self.task_id_counter, cmd_args))
        return Future(self.task_id_counter, cmd_args[-1], self.result_dict)
    
    def create_worker(self, gpu_id):
        return Worker(self.task_queue, self.task_status_dict, self.task_map, self.result_dict, gpu_id)

    def shutdown(self):
        logging.info(f"Shutting down TaskProcessingPool")
        self.stop_event.set()  # Signal to stop monitoring
        self.monitor_thread.join()
        with self.lock:
            for worker in self.workers:
                worker.stop_event.set()
            for worker in self.workers:
                worker.process.join(15)
                if worker.process.exitcode is None:
                    os.kill(worker.process.pid, signal.SIGKILL)  # Force kill the process
                    worker.process.join()
        logging.info(f"TaskProcessingPool shutdown complete")

class Future:
    def __init__(self, task_id, path, result_dict):
        self.task_id = task_id
        self.path = path
        self.result_dict = result_dict

    def result(self, timeout=1000):
        waiting_total = 0
        while self.task_id not in self.result_dict:
            if waiting_total > timeout:
                raise TimeoutError("Task {} timed out".format(self.task_id))
            waiting_total += 0.1
            time.sleep(0.1)  # Wait until the result is available
        return self.result_dict[self.task_id]
    
if __name__ == '__main__':
    import argparse
    import glob
    parser = argparse.ArgumentParser(description="Process FFT files with multiprocessing pool")
    parser.add_argument('-p', '--pool_size', type=int, required=True, help='Size of the pool')
    parser.add_argument('-d', '--directory', type=str, required=True, help='Directory containing FFT files')
    parser.add_argument('--zmax', type=str, required=True, help='zmax value')
    parser.add_argument('--wmax', type=str, required=True, help='wmax value')
    parser.add_argument('--sigma', type=str, required=True, help='sigma value')
    parser.add_argument('--numharm', type=str, required=True, help='numharm value')
    parser.add_argument('--batchsize', type=str, required=False, default='1', help='batchsize value')
    parser.add_argument('--pool_timeout', type=int, required=False, default=60*60, help='pool timeout in seconds')
    parser.add_argument('--gpu_count', type=int, required=False, default=1, help='Number of GPUs to use')
    # example:
    # python accelsearch_pipeline_cu.py -p 8 -d ffts --zmax 20 --wmax 0 --sigma 3.0 --numharm 16
    
    
    args = parser.parse_args()
    
    pool_size = args.pool_size
    directory = args.directory
    pool_timeout = args.pool_timeout
    fft_files = glob.glob(f'{directory}/*.fft')
    futures = []

    pool = Pool(pool_size, recreate_pool_timeout_seconds=pool_timeout, gpu_count=args.gpu_count)
    start = time.time()
    for fft_file in fft_files:
        cmd_args = [
            "accelsearch_pipeline_cu.py",  # Assumed script name
            "-zmax", args.zmax,
            "-wmax", args.wmax,
            "-sigma", args.sigma,
            "-numharm", args.numharm,
            "-batchsize", args.batchsize,
            fft_file
        ]
        futures.append(pool.submit(cmd_args))
    
    for future in futures:
        if future.result() != 0:
            logging.error(f"Task {future.path} failed with error code {future.result()}")
        else:
            logging.info(f"Task {future.path} completed successfully")
    time_elapsed = time.time() - start
    # logging
    logging.info(f"Time elapsed: {time_elapsed:.2f} seconds")
    logging.info(f"Processed {len(fft_files)} fft files")
    logging.info(f"Average time per fft file: {time_elapsed / len(fft_files):.2f} seconds")
    pool.shutdown()
    logging.info("All tasks completed")