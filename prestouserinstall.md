## 2. For PRESTO user
Because PrestoZL only add the cuda toolkit, all others are the same with PRESTO. If you have PRESTO(v4.0 preferred) running environment in your machine, you only need to install the CUDA toolkit and do some setting. 

**Attention：The PRESTO environment on your machine may not be compatible with PrestoZL**.

First Download and Install CUDA toolkit 11.8 From [here](https://developer.nvidia.com/cuda-11-8-0-download-archive).

Clone the repository to your local machine if you haven't already.
```
git clone https://github.com/zhejianglab/PrestoZL.git
```

Set the running path
```
vim ~/.bashrc
export PRESTO=/YourPathtoPrestoZL
export PATH=/YourPathtoPrestoZL/bin:${PATH}
export LD_LIBRARY_PATH=/YourPathtoPrestoZL/lib/:$LD_LIBRARY_PATH
source ~/.bashrc
```

Modify `Makefile` in `src`. 
```
CUDA_PATH       ?= /your/path/to/cuda
```
You can use default `GENCODE_FLAGS` setting, but if it goes wrong, set according to your GPU.

Compile PrestoZL
```
cd /YourPathtoPrestoZL/src
make makewisdom   % Only need to make when compiling for the first time
make cleaner && make
```

