#! /usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import subprocess
import sys

def readProfile(name):
    # Get file handle
    f = h5py.File(name, 'r')
    print(name)

    info = f['info']
    print('git commit: ', info['git-commit'][()][0])

    ranks = info['ranks'][()][0]
    print('ranks: ', ranks)

    timings = f['timings']
    collect = []
    for t in timings:
        count = timings[t+'/count'][()]
        time = timings[t+'/time'][()]
        # print(t, count, time)
        collect.append([t, ranks, count, time])

    return collect


if __name__ == '__main__':
    if(len(sys.argv) < 2):
        # Path
        path = '.'
    else:
        path = sys.argv[1]

    name = '/profile.hdf5'
    print(readProfile(path+name))
