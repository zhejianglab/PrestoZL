# PrestoZL
PrestoZL is a highly optimized, GPU-based pulsar search and analysis software developed by the team at the Astronomical Computing Research Center of Zhejiang Lab. It builds upon [Scott Ransom's PRESTO](https://github.com/scottransom/presto/tree/v4.0). The key difference between PrestoZL and Scott Ransom's PRESTO lies in the GPU optimization of the most time-consuming **"Jerk Search"** module, which has been tailored for GPU parallel processing pipelines.  **The search speed of PrestoZL can be accelerated by several tens of times faster than PRESTO, while maintaining the search logical and signal recovery capability fully consistent with PRESTO's Jerk Search.** 


<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure1.jpg" alt="Figure1" width="500">
  <p>Figure 1. Comparison of the Jerk Search Frameworks of PRESTO C and PrestoZL</p>
</div>

**Figure 1** compares the Jerk Search frameworks of PRESTO and PrestoZL. During each iteration of the r-step loop search, PrestoZL fuses the "ALL the harmonic summing and candidate search" logic into one GPU kernel, making the search process very efficient. At the same time, we have made a batch optimization to the r-step loop search. User can adjustment the parameters to achieve the maximum computational throughput according to your GPU.

**We also opensource a pipelined version of PrestZL，named PrestoZL-pipeline**, which eliminates the GPU stalls caused by extensive CPU computations. **Figure 2** illustrates the parallel framework of PrestoZL-pipeline. The framework enables a three-stage pipeline parallelism when processing consecutive FFT files within the same process. It effectively overlaps CPU computation time with GPU computation, and significantly improves the searching speed. The inner search logic is the PrestoZL.

<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure2.jpeg" alt="Figure2" width="600">
  <p>Figure 2. The three-stage pipeline framework of PrestoZL</p>
</div>

**We also opensource a GPU-accelerated version of De-dispersion** in `prepsubband_cu.c` with the performance far exceeds the CPU-based `prepsubband.c`. It can finish the entire de-dispersion process within half a minute, while keeping the results consistent with `prepsubband.c`. The logic and results are fully consistent with the code in [PRESTO's prepsubband.c](https://github.com/scottransom/presto/blob/v4.0/src/prepsubband.c).

**Figure 3** show the performance comparison results between different GPU implementations of PRESTO (a SOTA proprietary Presto GPU version, PrestoZL, and PrestoZL-pipeline) and PRESTO C. The metric used for comparison is the number of FFT files processed per minute on a single GPU. Both PrestoZL and PrestoZL-pipeline achieve significant performance improvements compare with PRESTO C.
<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure3.png" alt="Figure3" width="600">
  <p>Figure 3. Performance comparison between PRESTO C and different GPU versions under a single process</p>
</div>

*Test fits：FAST Dec+2737 02 05 arcdrift-M19 1443.fits,2GB; Test Environment: one A40 40G GPU，20 core CPU.

## For Pre PRESTO user(May Not Success)
If you have PRESTO(v4.0 prefered) running environment, all the environment dependencies are ready expect the CUDA toolkit. **Attention：The PRESTO environment set up on your machine may not be compatible with PrestoZL**.

First Download and Install CUDA toolkit V1.8 From [here](https://developer.nvidia.com/cuda-11-8-0-download-archive).

Compile PrestoZL
PrestoZL can be compiled as below:
```
cd /YourPathtoPrestoZL/src
make cleaner & make
```
Set the Runing Path
```
nano ~/.bashrc
export PATH="/YourPathtoPrestoZL/bin:${PATH}"
source ~/.bashrc
```

## Use Pre-built Docker Image(Recommand)

We have prepared a pre-built Docker image: **`zjlabastro/prestozl:latest`**，it is an environment-ready image to run PrestoZL. You can skip the `Build From Docker Image` steps below. Ensure you have Docker installed on your system. You can follow the instructions [here](https://github.com/zhejianglab/PrestoZL/blob/main/dockerInstall.MD) to install Docker. Once you have Docker installed on your system, the image can be fetched as follow:

Download pre-built docker image from dockerhub
```
docker pull zjlabastro/prestozl:latest
```
Run a Containter from the Image. -v can mount directories from host into a container, which can be used to share Fits data.
```
docker run -itd --name=prestozl_latest --gpus all --network=host -v /path/to/host/dir:/path/to/container/dir zjlabastro/prestozl:latest /bin/bash
```
Get into the Container and run PrestoZL.
```
docker exec -it prestozl_latest /bin/bash
```

## Build From Docker Image
You can simply build the PrestoZL enviroment from docker image. We have provided `Dockerfile`, and tested it on Ubuntu 20.04 with CUDA 11.7.1 .


### Prerequisites

Ensure you have Docker installed on your system. You can follow the instructions [here](https://github.com/zhejianglab/PrestoZL/blob/main/dockerInstall.MD) to install Docker.

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
make
```


## Usage
### PrestoZL 
`accelsearch_cu.c` serves as the entry point for the PrestoZL version of the jerk search program. The command has been expanded from the PRESTO C to include a batchsize parameter, which controls the number of rstep loops calculated on the GPU in each iteration. This parameter doesn't need to be explicitly set, its default value is **batchsize=8**. The lower the batchsize is ,the less GPU memory will be used. Other usages remain consistent with the PRRESTO . Here's an example:
```
accelsearch_cu -zmax 150 -wmax 150 -sigma 5.0 -numharm 16 -batchsize 2 tracking-M01_0047_DM9.15.fft
```
### PrestoZL-pipeline
To run PrestoZL-pipeline, you can use the python script at `bin/accelsearch_pipeline_cu.py`. `--pool_size` refers to the number of process to run concurrently, `--directory` refers to the input directory that stores the fft files, `--batchsize` is as the same meaning with PrestoZL. Here's an example:
```
accelsearch_pipeline_cu.py --pool_size 8 --directory ffts --zmax 150 --wmax 150 --sigma 3.0 --numharm 16 --batchsize 2
```
### De-dispersion
To run the GPU-accelerated version of de-dispersion, you can use the command `prepsubband_cu`, other arguments are the same with the `prepsubband` used in PRESTO C. Here's an example:
```
prepsubband_cu -nobary -numout 262144 -nsub 3280 -lodm 59.4 -dmstep 0.1 -numdms 564 -downsamp 1 -mask test1_rfifind.mask -o ./psb/FRB121102_tracking-M01_0706_ds1_0 FRB121102_tracking-M01_0706_ds1_0.fits
```
The method of running other parts of PRESTO is the same with PRESTO C.

## Acknowledgement
The following individuals contributed to this project, (listed in alphabetical order): Chen Huaxi, Yinan Ke, Mao Kuang, Pan Qiuhong, Tang Zhaorong, Tang Xuefei, Wang Qi, Wang Pei

If you have any question, be free to contact me:  maok@zhejianglab.org

If you use the code in your experiment, please cite this repo, thank you.

```
@misc{PrestoZL,
  author = {Zhejiang Lab},
  title = {PrestoZL，GPU version of PRESTO},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/zhejianglab/PrestoZL.git}},
}
```
