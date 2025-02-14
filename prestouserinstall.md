## 2. For PRESTO user
There are two solutions. 

### 2.1 First solution:Start with PrestoZL source code
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

### 2.2 Second solution: copy files from FrestoZL to your PRESTO source code
PrestoZL is based on Scott Ransom's PRESTO v4.0. If your existing PRESTO environment is a newer version and you wish to continue using it, You can simply copy the following files from PrestoZL into your PRESTO C's corresponding directory and then recompile. The required changes are as follows:

1. Add new files from PrestoZL to your PRESTO source.
```
bin/accelsearch_pipeline_cu.py
include/cuda_helper.h
include/exception.h
include/helper_cuda.h
include/helper_cuda_drvapi.h
include/helper_cusolver.h
include/helper_functions.h
include/helper_gl.h
include/helper_image.h
include/helper_math.h
include/helper_multiprocess.h
include/helper_nvJPEG.hxx
include/helper_string.h
include/helper_timer.h
src/accel_kernels.cu
src/accel_utils_cu.c
src/accelsearch_cu.c
src/accelsearch_lib.c
src/accelsearch_maxsize.c
src/prepsubband_cu.c
src/prepsubband_kernels.cu
src/prepsubband_maxsize.c
```

2. Overwrite existing files in your PRESTO source with the corresponding files from PrestoZL.
```
include/accel.h
include/accelsearch_cmd.h
src/accel_utils.c
src/accelsearch_cmd.c
src/characteristics.c
src/dcdflib.c
src/fftcalls.c
```
3. Download and Install CUDA toolkit 11.8 From [here](https://developer.nvidia.com/cuda-11-8-0-download-archive).

4. Modify existing files according to the changes in PrestoZL.
```
src/Makefile
```
Update your src/Makefile to include the CUDA environment configuration and the compile configuration for `accelsearch_cu` and `prepsubband_cu` according to PrestoZL's src/Makefile.

5. Run `make` again in your PRESTO/src directory.
```
cd /YourPathtoPrestoZL/src
make cleaner && make
```
