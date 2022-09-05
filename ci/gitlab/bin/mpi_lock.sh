#!/bin/bash

#
# Utility script for serial commands when using multiple processes
#

# get a hash based on args to create a unique lock
HASHMD5=`echo "$@" | md5sum | cut -f1 -d" "`
lockfile=$CI_CACHE_FOLDER/.lock-mpi.$HASHMD5
tmpfile=${lockfile}.$$
echo $$ > $tmpfile
# ln is atomic
if ln $tmpfile $lockfile 2>&-; then
    # unlock and cleanup on exit
    trap 'rm -f $tmpfile $lockfile' 0
    echo locked
    eval "$@"
    # forward return value
    retVal=$?
    echo unlocked
    if [ $retVal -ne 0 ]; then
        exit $retVal
    fi
else
    # cleanup on exit
    trap 'rm -f $tmpfile' 0
    echo locked by $(<$lockfile)
    while [ -f $lockfile ]; do
        # to not hammer file system
        sleep 0.1s
    done
    echo $$ resuming
fi
