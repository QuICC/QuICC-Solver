FROM finkandreas/ubuntu-quicc-buildbase:20.04

COPY . /QuICC.src

RUN mkdir -p /QuICC.src/build \
  && cd /QuICC.src/build \
  && cmake .. -DCMAKE_BUILD_TYPE=RelWithAssert -DCMAKE_MESSAGE_LOG_LEVEL=VERBOSE \
    -DQUICC_TESTSUITE_POLYNOMIAL=ON \
    -DQUICC_TESTSUITE_FRAMEWORK=ON \
    -DQUICC_TESTSUITE_TRANSFORM=ON \
    -DQUICC_TESTSUITE_PROFILING=ON \
    -DQUICC_TESTSUITE_SPARSESM=ON \
  && make -j $(grep processor /proc/cpuinfo | wc -l)
