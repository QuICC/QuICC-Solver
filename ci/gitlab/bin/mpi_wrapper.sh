#!/bin/bash

#
# Utility script to wrap non executables for mpi launch
#


eval "$@"
retVal=$?
if [ $retVal -ne 0 ]; then
    exit $retVal
fi
