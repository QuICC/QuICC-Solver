# Use docker image

In order to simplify setting up a local version of QuICC, a docker image with required dependencies installed is provided. The latest image can be found on docker hub at [quicc/buildbase](https://hub.docker.com/repository/docker/quicc/buildbase/general). The setup described below will mount the `</path/to/QuICC_mnt>` directory into the container. Everything that needs to remain even after stopping the docker container should be stored in this directory. For example, source and build directory should be created in `</path/to/QuICC_mnt>`. In the commands below, `</path/to/QuICC_mnt>` needs to be replaced with an actual path.

## Linux and MacOS

### Setup docker

The following steps should produce a working environment:
 1. Install docker:
    - *Linux:* Follow instructions for your distribution.
    - *MacOS:* Install docker desktop from [Install on Mac](https://docs.docker.com/desktop/install/mac-install/).
 2. Add your SSH keys to the SSH agent if not already present, i.e. for each key you need, run:
    ```bash
    ssh-add </path/to/SSH_KEY>
    ```
 3. Start docker
 4. Open Terminal:
    - Pull docker image from hub
     ```bash
     docker pull quicc/buildbase:latest
     ```
    - Create directory where sources and builds will be kept. This directory is persistent and will be mounted into the container.
    ```bash
    mkdir </path/to/QuICC_mnt>
    ```
    - Clone QuICC repository into `<path/to/QuICC_mnt>`, for example to clone the `main` branch:
    ```bash
    git clone -b main git@github.com:QuICC/QuICC </path/to/QuICC_mnt>/QuICC
    ```
    
### Run docker image in Linux

The directory `</path/to/QuICC_mnt>` will be mounted into the container under `/QuICC`. The content of this directory is persistent and will not disappear when the container is deleted.
```bash
docker run --rm -it \
-v </path/to/QuICC_mnt>:/QuICC \
--mount type=bind,src="$SSH_AUTH_SOCK",target=/run/host-services/ssh-auth.sock \
-e SSH_AUTH_SOCK="/run/host-services/ssh-auth.sock" \
quicc/buildbase:latest
```

If everything worked, you should now be inside the container. Create a build directory in `/QuICC` and proceed as usual to configure, build and run QuICC.
    
### Run docker image in MacOS

The directory `</path/to/QuICC_mnt>` will be mounted into the container under `/QuICC`. The content of this directory is persistent and will not disappear when the container is deleted.
```bash
docker run --rm -it \
-v </path/to/QuICC_mnt>:/QuICC \
--mount type=bind,src="/run/host-services/ssh-auth.sock",target=/run/host-services/ssh-auth.sock \
-e SSH_AUTH_SOCK="/run/host-services/ssh-auth.sock" \
quicc/buildbase:latest
```

If everything worked, you should now be inside the container. Create a build directory in `/QuICC` and proceed as usual to configure, build and run QuICC.

# How to run on Piz-Daint

## Cpu
```sh
module load daint-mc
module load sarus
export DOCKERHUB_USERNAME=<your-user-name>
srun -A s1111 -C mc sarus pull $DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:22.04
srun -A s1111 -C mc --pty sarus run -t --mount=type=bind,source=$SCRATCH,destination=$SCRATCH $DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:22.04 bash
```

## Gpu
```sh
module load daint-gpu
module load sarus
export DOCKERHUB_USERNAME=<your-user-name>
srun -A c19 -C gpu sarus pull $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04
srun -A c19 -C gpu --pty sarus run -t --mount=type=bind,source=$SCRATCH,destination=$SCRATCH $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04 bash
```

# How to build base images

## Cpu spack baseimage
```sh
export DOCKERHUB_USERNAME=<your-user-name>
docker build --pull --force-rm -f <path-to-quicc>/ci/docker/baseimage/Dockerfile_spack_baseimage_cpu -t $DOCKERHUB_USERNAME/ubuntu-spack:22.04-0.18.0 .
docker push $DOCKERHUB_USERNAME/ubuntu-spack:22.04-0.18.0
```
## Cpu QuICC baseimage
```sh
docker build --pull --force-rm --build-arg CSCS_REGISTRY_PATH=$DOCKERHUB_USERNAME --build-arg NUM_PROCS=2 --build-arg TARGET=none -f <path-to-quicc>/ci/docker/baseimage/Dockerfile_quicc_baseimage_cpu -t $DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:22.04 .
docker push $DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:22.04
```
## Cpu QuICC linear stability
```sh
docker build --pull --force-rm --build-arg CSCS_REGISTRY_PATH=$DOCKERHUB_USERNAME --build-arg NUM_PROCS=2 --build-arg TARGET=none -f <path-to-quicc>/ci/docker/baseimage/Dockerfile_quicc_stability_cpu -t $DOCKERHUB_USERNAME/ubuntu-quicc-stability:22.04 .
docker push $DOCKERHUB_USERNAME/ubuntu-quicc-stability:22.04
```

## Gpu spack baseimage
```sh
export DOCKERHUB_USERNAME=<your-user-name>
docker build --pull --force-rm -f <path-to-quicc>/ci/docker/baseimage/Dockerfile_spack_baseimage_gpu -t $DOCKERHUB_USERNAME/cuda-ubuntu-spack:11.7.1-22.04-0.18.0 .
docker push $DOCKERHUB_USERNAME/cuda-ubuntu-spack:11.7.1-22.04-0.18.0
```
## Gpu QuICC baseimage
```sh
docker build --pull --force-rm --build-arg CSCS_REGISTRY_PATH=$DOCKERHUB_USERNAME --build-arg NUM_PROCS=2 --build-arg TARGET=none -f <path-to-quicc>/ci/docker/baseimage/Dockerfile_quicc_baseimage_gpu -t $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04 .
docker push $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04
```

# How to make multi-architecture images

Tag and push the images:

```sh
export DOCKERHUB_USERNAME=<your-user-name>
docker push $DOCKERHUB_USERNAME/buildbase:amd64-broadwell
docker push $DOCKERHUB_USERNAME/buildbase:arm64-v8
```

Create the new manifest:

```sh
docker manifest create \
$DOCKERHUB_USERNAME/buildbase:latest \
--amend $DOCKERHUB_USERNAME/buildbase:amd64-broadwell \
--amend $DOCKERHUB_USERNAME/buildbase:arm64-v8
```

Push the new manifest:

```sh
docker manifest push $DOCKERHUB_USERNAME/buildbase:latest
```
