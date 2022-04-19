# Build docker image for QuICC

```bash
cd /path/to/QuICC
docker build --pull --force-rm -f docker/baseimage/Dockerfile_quicc_baseimage -t DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:20.04 .
docker push DOCKERHUB_USERNAME/ubuntu-quicc-buildbase:20.04
```

# Navigate baseimage container

```bash
docker pull finkandreas/ubuntu-quicc-buildbase:20.04
docker run --rm -it finkandreas/ubuntu-quicc-buildbase:20.04
```
