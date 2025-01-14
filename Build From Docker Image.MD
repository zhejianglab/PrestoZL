## 3. Build From Docker Image
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
