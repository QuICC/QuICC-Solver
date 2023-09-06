#!/bin/bash

#
# Utility script to wrap reframe for mpi launch
# only rank 1 calls reframe
#

if [ -n "$OMPI_COMM_WORLD_RANK" ]; then
    if [ $OMPI_COMM_WORLD_RANK -eq 0 ]; then
        echo "OPEN-MPI detected"
    fi
    RANK=$OMPI_COMM_WORLD_RANK
elif [ -n "$SLURM_PROCID" ] ; then
    if [ $SLURM_PROCID -eq 0 ]; then
        echo "SLURM detected"
    fi
    RANK=$SLURM_PROCID
else
    echo "no known mpi enviroment detected"
    exit 1
fi

if [ $RANK -eq 0 ]; then
    eval "$1"
else
    eval "$2"
fi

retVal=$?
if [ $retVal -ne 0 ]; then
    exit $retVal
fi
