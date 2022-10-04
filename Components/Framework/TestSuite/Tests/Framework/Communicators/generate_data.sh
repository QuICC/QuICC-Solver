#!/bin/bash

d1d=9
d2d=19
d3d=19

schemes=("WLFl" "WLFm" "SLFl" "SLFm")
nps=("4")
algorithms=("tubular" "single1d" "single2d")

# Give test executable as first argument to script

mpinp="mpirun -np "
exe=$(realpath $1)
for alg in ${algorithms[@]};
do
  mkdir -p ${alg};
  pushd ${alg};
  for sch in ${schemes[@]};
  do
    mkdir -p ${sch};
    pushd ${sch};
    for np in ${nps[@]};
    do
      mkdir -p Np${np};
      pushd Np${np};
      rm -f *.gxl *.vtp
      ${mpinp} ${np} ${exe} [${sch}] --algorithm ${alg} --dim1D ${d1d} --dim2D ${d2d} --dim3D ${d3d}
      popd;
    done
    popd;
  done
  popd
done
