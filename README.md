# PrestoZL
PrestoZL is a highly optimized GPU-based pulsar search and analysis software developed by the team at the Astronomical Computing Research Center of Zhejiang Lab. It is developed based on [Scott Ransom's Presto C code](https://github.com/scottransom/presto). The key difference between PrestoZL and Scott Ransom's Presto C is the GPU optimization of the most time-consuming "Jerk Search" module, which is tailored for GPU parallel processing pipelines. This optimization significantly accelerates the search process. As the search parameters `(zmax, wmax)` increase, the performance improvement becomes more noticeable. **PrestoZL's search results is equivalent to the Scott Ransom's PRESTO C.**

<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure1.jpg" alt="Figure1" width="500">
  <p>Figure 1. Comparison of the Jerk Search Frameworks of PRESTO C and PrestoZL</p>
</div>

**Figure 1** compares the jerk search frameworks of PRESTO C and PrestoZL. During each iteration of the r-step loop search, PrestoZL combines the "harmonic summing and candidate search" logic. This combination allows us to complete each round of search within one cuda kernel, making the search process very efficient. At the same time, we have made a batch modification to the r-step loop search. This means that, depending on the GPU's memory capacity and different search parameters, adjustments can be made to achieve the maximum computational throughput.

**We also provide a pipelined version of PrestZL，named PrestoZL-pipeline**, which eliminates the GPU stalls caused by extensive CPU computations. **Figure 2** illustrates the parallel framework of PrestoZL-pipeline. The framework enables a three-stage pipeline parallelism when processing consecutive FFT files within the same process. It effectively overlaps CPU computation time with GPU computation, and significantly improving the searching speed. The inner search logic is the PrestoZL.

<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure2.jpeg" alt="Figure2" width="600">
  <p>Figure 2. The three-stage pipeline framework of PrestoZL</p>
</div>

**We also provide a GPU-accelerated version of De-dispersion** in `prepsubband_cu.c` with the performance far exceeds the CPU-based `prepsubband.c`. It can finish the entire de-dispersion process within half a minute, while keeping the results consistent with `prepsubband.c`.

**Figure 3** and **Figure 4** show the performance comparison results between different GPU implementations of PRESTO (a SOTA proprietary Presto GPU version, PrestoZL, and the pipelined version of PrestoZL) and PRESTO C. The metric used for comparison is the number of FFT files processed per minute under a single process and eight concurrent processes on a single GPU. Optimizing GPU programs for multi-process search on the same GPU is challenging, requiring fine-tuned optimization and utilization of GPU computational resources and memory access. Both PrestoZL and the pipelined version of PrestoZL achieve significant performance improvements compare with PRESTO C.

<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure3.png" alt="Figure3" width="600">
  <p>Figure 3. Performance comparison between PRESTO C and different GPU versions under a single process</p>
</div>

<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure4.png" alt="Figure4" width="600">
  <p>Figure 4. Performance comparison between PRESTO C and different GPU versions under 8 concurrent processes</p>
</div>


## Build From Docker Image
You can simply build the PrestoZL enviroment from docker image. We have provided `Dockerfile`, and tested it on Ubuntu 20.04 with CUDA 11.7.1 .

### Prerequisites

Ensure you have Docker installed on your system. You can follow the instructions [here](https://docs.docker.com/get-docker/) to install Docker.

### Steps to Build the Docker Image

1. **Clone the Repository**:

   First, clone the repository to your local machine if you haven't already.
   ```
   git clone https://github.com/zhejianglab/PrestoZL.git
   cd PrestoZL
   ```
2. **Build the Docker Image**:

   Use the following command to build the Docker image. Make sure you are in the directory where the Dockerfile is located.
   ```
   docker build -t prestozl:latest .
   ```
   This command will build the Docker image and tag it as `prestozl:latest`.
3. **Run the Docker Container**:

   Once the image is built, you can run a container from it. Here is an example command to run the container interactively.
   ```
   docker run -itd --name=prestozl_latest --gpus all --network=host prestozl:latest /bin/bash
   docker exec -it prestozl_latest /bin/bash
   ```
   These commands will create a Docker container named `prestozl_latest` and enter it.

### Compile PrestoZL
PrestoZL can be compiled in the container as below commands:
```
cd /home/soft/presto/src
make makewisdom
make
```
The command `make makewisdom` only needs to be executed during the first compilation.


## Usage
### PrestoZL 
`accelsearch_cu.c` serves as the entry point for the PrestoZL version of the jerk search program. The command has been expanded from the Presto C to include a batchsize parameter, which controls the number of rstep loops calculated on the GPU in each iteration. This parameter doesn't need to be explicitly set, its default value is **batchsize=8**. The lower the batchsize is ,the less GPU memory will be used. Other usages remain consistent with the Presto C . Here's an example:
```
accelsearch_cu -zmax 200 -wmax 500 -sigma 5.0 -numharm 16 -batchsize 2 tracking-M01_0047_DM9.15.fft
```
### PrestoZL-pipeline
To run PrestoZL-pipeline, you can run a python script at `bin/accelsearch_pipeline_cu.py`. `--pool_size` refers to the number of process to run concurrently, `--directory` refers to the input directory that stores the fft files, `--batchsize` is as the same meaning with PrestoZL. Here's an example:
```
python accelsearch_pipeline_cu.py --pool_size 8 --directory ffts --zmax 20 --wmax 0 --sigma 3.0 --numharm 16 --batchsize 2
```
### De-dispersion
To run the GPU-accelerated version of de-dispersion, you can use the command `prepsubband_cu`, other arguments are the same as the `prepsubband` used in PRESTO C. Here's an example:
```
prepsubband_cu -nobary -numout 262144 -nsub 3280 -lodm 59.4 -dmstep 0.1 -numdms 564 -downsamp 1 -mask test1_rfifind.mask -o ./psb/FRB121102_tracking-M01_0706_ds1_0 FRB121102_tracking-M01_0706_ds1_0.fits
```
The method of running other parts of PRESTO is the same with PRESTO C.

## Acknowledgement
The following individuals contributed to this project, (listed in alphabetical order): Chen Huaxi, Yinan Ke, Mao Kuang, Pan Qiuhong, Tang Zhaorong, Tang Xuefei, Wang Qi, Wang Pei

If you have any question, be free to contact me:  maok@zhejianglab.com

If you use the code in your experiment, please cite this repo, thank you.
```
@misc{PrestoZL,
 author = {Zhejiang Lab},
 title = {PrestoZL，Zhejiang Lab Presto GPU},
 year = {2024},
 howpublished = {\url{https://github.com/zhejianglab/PrestoZL.git}},
 note = {Accessed: 2024-06-24}
}
```
