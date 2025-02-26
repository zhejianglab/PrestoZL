# PrestoZL
PrestoZL is a highly optimized, GPU-based pulsar search and analysis software developed by the team at the Astronomical Computing Research Center of Zhejiang Lab. It was developped based on the Scott Ransom's PRESTO v4.0. The key difference between them lies in the GPU optimization of the most time-consuming **"Jerk Search"** module, which has been tailored for GPU parallel processing pipelines.  **The search speed of the GPU-based Jerk Search in PrestoZL can be accelerated by several tens of times faster than CPU-based Jerk Search in PRESTO, while maintaining the search efficacy and signal recovery capability fully consistent with [CPU-based Jerk Search in PRESTO v4.0](https://github.com/scottransom/presto/tree/v4.0).** 

PrestoZL is particularly suitable for searching pulsars using large (zmax, wmax) parameters with long observation times. Currently, we have tested our code on observation files of approximately one hour with zmax=200 and wmax=500. If you encounter any issues during the processing of larger parameters or longer observation times, please contact us.

<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure1.jpg" alt="Figure1" width="500">
  <p>Figure 1. Comparison of the Jerk Search Frameworks of PRESTO C and PrestoZL</p>
</div>

**Figure 1** compares the Jerk Search frameworks of PRESTO C and PrestoZL. During each iteration of the `r-step` loop, PrestoZL fuses the "ALL the harmonic summing and candidate search" logic into one GPU kernel, making the search process very efficient. We also support batch the calculation of several `rstep` into one. User can adjustment the `batchsize` parameters to achieve the maximum computational throughput according to your GPU.

**We also opensource a pipelined version of PrestZL，named PrestoZL-pipeline**, which eliminates the GPU stalls caused by extensive CPU computations. **Figure 2** illustrates the parallel framework of PrestoZL-pipeline. The framework enables a three-stage pipeline parallelism when processing **consecutive FFT files** within the same process. It effectively overlaps CPU computation time with GPU computation among next, current and previous processing FFT files. The inner search logic is the PrestoZL.

<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure2.jpeg" alt="Figure2" width="600">
  <p>Figure 2. The three-stage pipeline framework of PrestoZL</p>
</div>

**We also opensource a GPU-accelerated version of De-dispersion** in `prepsubband_cu.c` with the performance far exceeds the CPU-based `prepsubband.c`. It can be 100 times faster than the CPU based de-dispersion version, and the processing results are fully consistent with the code in [PRESTO's prepsubband.c](https://github.com/scottransom/presto/blob/v4.0/src/prepsubband.c).

**Figure 3** show the performance comparison results between different GPU implementations of PRESTO (a SOTA proprietary Presto GPU version, PrestoZL, and PrestoZL-pipeline) and PRESTO C. The metric used for comparison is the number of FFT files processed per minute on a single GPU. Both PrestoZL and PrestoZL-pipeline achieve significant performance improvements compare with PRESTO C.
<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure3.png" alt="Figure3" width="600">
  <p>Figure 3. Performance comparison between PRESTO C and different GPU versions under a single process</p>
</div>

*Test fits：FAST Dec+2737 02 05 arcdrift-M19 1443.fits,2GB; Test Environment: one A40 40G GPU，20 core CPU.

## PrestoZL Environment Setting
There are three ways to set up the PrestoZL running environment. 
### 1. Use Pre-built Docker Image(Recommand)

We have prepared a pre-built Docker image: **`zjlabastro/prestozl:v1.0`**，it is an environment-ready image to run PrestoZL. Ensure you have Docker and nvidia-container-toolkit installed on your system, you can follow the instructions [here](https://github.com/zhejianglab/PrestoZL/blob/main/dockerInstall.MD). Then the image can be fetched as follow:

Download pre-built docker image from dockerhub. 
```
docker pull zjlabastro/prestozl:v1.0
```
*If you have trouble in visiting dockerhub, download the Image from [BaiduNetdisk](https://pan.baidu.com/s/11EMfbEFDQdRxn-tpA8EWrw?pwd=jhd3), and follow readme.txt.

Run a Containter from the Image. `-v` mounts directories from host file system into the docker container for data share.
```
docker run -itd --name=prestozl_latest --gpus all --network=host -v /path/to/your/host/dir:/home/soft/prestodata zjlabastro/prestozl:v1.0 /bin/bash
```
Get into the Container and run PrestoZL.
```
docker exec -it prestozl_latest /bin/bash
```

### 2. For PRESTO C user without Docker
If you have PRESTO v4.0 environment on your machine, and really DO NOT want to use Docker. Try the solution [here](https://github.com/zhejianglab/PrestoZL/blob/main/prestouserinstall.md).


### 3. Build From Dockerfile
We have provided `Dockerfile` to build the PrestoZL Image by yourself. Follow the instruction [here](https://github.com/zhejianglab/PrestoZL/blob/main/Build%20From%20Docker%20Image.MD).


## Usage
### PrestoZL 
`accelsearch_cu.c` serves as the entry point for the PrestoZL's Jerk Search. The command has been expanded from the PRESTO to include a `batchsize` parameter, which batching the calculation of several `rstep` in each while iteration. This parameter doesn't need to be explicitly set, its default value is **batchsize=8**. The lower the batchsize is ,the less GPU memory will be used. Other usages remain consistent with the PRRESTO . Here's an example:
```
accelsearch_cu -zmax 50 -wmax 50 -sigma 5.0 -numharm 16 -batchsize 8 yourFFTfile.fft
```
To use multiple processes(set by `-P`) to process many FFT files in `/yourPathtoFFTfiles` concurrently:
```
ls /yourPathtoFFTfiles/*.fft |  xargs -P 8 -n 1 accelsearch_cu -zmax 50 -wmax 50 -sigma 3.0 -numharm 16 -batchsize 8
```
### PrestoZL-pipeline
To run PrestoZL-pipeline, you can use the python script at `bin/accelsearch_pipeline_cu.py`. `--pool_size` refers to the number of concurrent running process in a GPU, each process is an FFT file processing pipeline, `--directory` refers to the directory of the input fft files, `--batchsize` is as the same meaning with PrestoZL. Here's an example:
```
accelsearch_pipeline_cu.py --pool_size 8 --directory /yourPathtoFFTfiles --zmax 50 --wmax 50 --sigma 3.0 --numharm 16 --batchsize 8
```
### De-dispersion
To run the GPU-accelerated version of de-dispersion, you can use the command `prepsubband_cu`, other arguments are the same with the `prepsubband` used in PRESTO. Here's an example:
```
prepsubband_cu -nobary -numout 262144 -nsub 3280 -lodm 59.4 -dmstep 0.1 -numdms 564 -downsamp 1 -mask your_rfifind.mask -o /path/to/outputdir/withprefix /yourpathtoFITSfiles/my.fits
```
The usage of other commands are the same with PRESTO.

## Acknowledgement
The following individuals contributed to this project, (listed in alphabetical order): Chen Huaxi, Mao Kuang, Pan Qiuhong, Tang Zhaorong, Tang Xuefei, Wang Qi, Wang Pei

If you have any question, be free to contact me:  maok@zhejianglab.org

If you use the code in your experiment, please cite this repo, thank you.

```
@misc{PrestoZL,
  author = {ZhejiangLab},
  title = {PrestoZL GPU version of PRESTO},
  publisher = {GitHub},
  url = {https://github.com/zhejianglab/PrestoZL.git},
  year = {2024}
}
```
