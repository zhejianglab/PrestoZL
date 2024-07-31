# PrestoZL
PrestoZL is a highly optimized GPU-based pulsar search and analysis software developed by the team at the Astronomical Computing Research Center of Zhejiang Lab. It is developed based on [Scott Ransom's Presto C code](https://github.com/scottransom/presto). The difference between this software and Scott Ransom's Presto C lies in the gpu optimization of the most time-consuming Fourier-domain acceleration search module to suit GPU parallel processing pipelines, which significantly accelerating the search procedual. The larger the search parameters `(zmax, wmax)`, the more noticeable the performance improvement, and **the search result is equivalent to the Scott Ransom's Presto C**.

In **Figure 1**, the GPU version of Fourier-domain acceleration（jerk search） is compared to the logic of accelsearch.c in Presto C. During each iteration of the r-step loop search, we combine the "harmonic summing and candidate search" logic. This combination allows us to complete each round of search within one Cuda Kernel, making the search process very efficient. At the same time, we have made a batch modification to the r-step loop search. This means that, depending on the GPU's memory capacity and different search parameters, adjustments can be made to achieve the maximum computational throughput.
<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure1.jpg" alt="Figure1" width="600">
  <p>Figure 1. Comparison of the program frameworks of Presto C and PrestoZL in the accelsearch stage</p>
</div>

**Figure 2** shows the performance comparison between our GPU accelerated search algorithm, Presto C, and a SOTA proprietary Presto GPU version. PrestoZL demonstrates optimal processing throughput. In the test using a single A40 GPU for accelerated search, because of the well-optimized memory usage, our version supports accelerated search run with 12 concurrent processes on one GPU(A40 40G) at search parameters of wmax=150 and zmax=150.
<div align="center">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure2.png" alt="Figure2" width="600">
  <img src="https://github.com/zhejianglab/PrestoZL/raw/main/resource/Figure3.png" alt="Figure3" width="600">
  <p>Figure 2. Performance comparison under single process and multi-process</p>
</div>

## Build From Docker Image
You can simply build the PrestoZJ enviroment from docker image. We have provided `Dockerfile`, and tested it on Ubuntu 20.04 with CUDA 11.7.1 .

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
   docker build -t prestozj:latest .
   ```
   This command will build the Docker image and tag it as `prestozj:latest`.
3. **Run the Docker Container**:

   Once the image is built, you can run a container from it. Here is an example command to run the container interactively.
   ```
   docker run -itd --name=prestozj_latest --gpus all --network=host prestozj:latest /bin/bash
   docker exec -it prestozj_latest /bin/bash
   ```
   These commands will create a Docker container named `prestozj_latest` and enter it.


## Usage
The Fourier-domain acceleration (accelsearch_cu.c) section serves as the entry point for the PrestoZL version of the Fourier-domain acceleration program. The command has been expanded from the Presto C to include a batchsize parameter, which controls the number of rstep loops calculated on the GPU in each iteration. This parameter doesn't need to be explicitly set, its default value is **batchsize=8**. The lower the batchsize is ,the less GPU memory will be used. Other usages remain consistent with the Presto C . Here's an example:
```
accelsearch_cu -zmax 200 -wmax 500 -sigma 5.0 -numharm 16 -batchsize 2 tracking-M01_0047_DM9.15.fft
```
The method of running other parts of Presto is the same with Presto C.

## Acknowledgement
The following individuals contributed to this project, (listed in alphabetical order): Chen Huaxi, Ke Yinan, Mao Kuang, Pan Qiuhong, Tang Zhaorong, Tang Xuefei, Wang Qi, Wang Pei

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
