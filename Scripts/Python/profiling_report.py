#! /usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import subprocess
import sys
import numpy as np

def reportTimings(prefix, timings, max_lvl, lvl = 1):
    top_timings = list()
    if prefix not in timings and prefix:
        print('\t'*(lvl-1) + prefix.split('-')[-1])
    for t in timings:
        if t.startswith(prefix + '-'):
            lt = t.split('-')
            nt = len(lt)
            if nt > lvl and nt <= max_lvl:
                tt = lt[0:lvl+1]
                if len(tt) > 1:
                    tt = '-'.join(tt)
                else:
                    tt = tt[0]
                if tt not in top_timings:
                    top_timings.append(tt)
        elif t == prefix:
            formatTiming(t, timings[t], lvl-1)
    for t in top_timings:
        reportTimings(t, timings, max_lvl, lvl + 1)

def formatTiming(t, timings, tabs):
    count = timings['count'][()]
    ts = timings['time'][()]
    t_min = np.min(ts)
    t_avg = np.average(ts)
    t_max = np.max(ts)
    opName = t.split('-')[-1]
    print('\t'*tabs+f'{opName+":":<30} {t_min:.2e} / {t_avg:.2e} / {t_max:.2e}')

def readProfile(name):
    # Get file handle
    f = h5py.File(name, 'r')
    print(name)

    info = f['info']
    print('git commit: ', info['git-commit'][()][0])

    ranks = info['ranks'][()][0]
    print('ranks: ', ranks)

    timings = f['timings']

    print('Timing format: min / avg / max')
    # Forward
    reportTimings('Fwd', timings, 3)

    # Backward
    reportTimings('Bwd', timings, 3)

    # Trivial
    reportTimings('Trivial', timings, 2)

    # Diagnostic
    reportTimings('Diagnostic', timings, 2)

    # Prognostic
    reportTimings('Prognostic', timings, 2)

    # Timestep
    reportTimings('Timestep', timings, 2)

    # Control
    reportTimings('Control', timings, 2)

    # IO
    reportTimings('IO', timings, 2)

    # Worland operators
    reportTimings('Worland::Poly', timings, 3)
    reportTimings('Worland::Fft', timings, 3)


if __name__ == '__main__':
    if(len(sys.argv) < 2):
        # Path
        path = '.'
    else:
        path = sys.argv[1]

    name = '/profile.hdf5'
    readProfile(path+name)
