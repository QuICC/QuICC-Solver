#!/bin/bash

patch=$1

exec 3>&1                    # Save the place that stdout (1) points to.
out=$(git apply --check $patch 2>&1 1>&3)  # Run command.  stderr is captured.
exec 3>&-                    # Close FD #3.

if [[ "$out" =~ .*"error".* ]]  ; then
    echo Patch is going to fail, skip
else
    echo Patch can be safely applied
    exec git apply $patch
fi
