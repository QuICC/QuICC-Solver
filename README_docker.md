# How to run on Piz-Daint

## Cpu
```sh
module load daint-mc
module load sarus
export DOCKERHUB_USERNAME=giacasti
srun -A s1111 -C mc sarus pull $DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:22.04
srun -A s1111 -C mc --pty sarus run -t --mount=type=bind,source=$SCRATCH,destination=$SCRATCH $DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:22.04 bash
```

## Gpu
```sh
module load daint-gpu
module load sarus
export DOCKERHUB_USERNAME=giacasti
srun -A c19 -C gpu sarus pull $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04
srun -A c19 -C gpu --pty sarus run -t --mount=type=bind,source=$SCRATCH,destination=$SCRATCH $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04 bash
```

# How to build base images

## Cpu spack baseimage
```sh
export DOCKERHUB_USERNAME=<your-user-name>
docker build --pull --force-rm -f ci/docker/baseimage/Dockerfile_spack_baseimage_cpu -t $DOCKERHUB_USERNAME/ubuntu-spack:22.04-0.18.0 .
docker push $DOCKERHUB_USERNAME/ubuntu-spack:22.04-0.18.0
```
## Cpu QuICC baseimage
```sh
docker build --pull --force-rm --build-arg CSCS_REGISTRY_PATH=$DOCKERHUB_USERNAME --build-arg NUM_PROCS=2 -f ci/docker/baseimage/Dockerfile_quicc_baseimage_cpu -t $DOCKERHUB_USERNAME/ubuntu-spack:22.04-0.18.0 .
docker push $DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:22.04
```

## Gpu spack baseimage
```sh
export DOCKERHUB_USERNAME=<your-user-name>
docker build --pull --force-rm -f ci/docker/baseimage/Dockerfile_spack_baseimage_gpu -t $DOCKERHUB_USERNAME/cuda-ubuntu-spack:11.7.1-22.04-0.18.0 .
docker push $DOCKERHUB_USERNAME/cuda-ubuntu-spack:11.7.1-22.04-0.18.0
```
## Cpu QuICC baseimage
```sh
docker build --pull --force-rm --build-arg CSCS_REGISTRY_PATH=$DOCKERHUB_USERNAME --build-arg NUM_PROCS=2 -f ci/docker/baseimage/Dockerfile_quicc_baseimage_gpu -t $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04 .
docker push $DOCKERHUB_USERNAME/cuda-ubuntu-quicc-buildbase:11.7.1-22.04
```

