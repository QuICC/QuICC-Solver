#! /usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import subprocess
import sys
import numpy as np

def printHeader(name, tabs):
    print('\t'*tabs+f'{name+":":<20}')

def printTiming(name, tag, db, tabs, debug = True, global_tot = None):
    t_tot = None
    if tag in db:
        timings = db[tag]
        count = timings['count'][()]
        ts = timings['time'][()]
        t_tot = np.average(ts)*count
        t_min = np.min(ts)
        t_avg = np.average(ts)
        t_max = np.max(ts)
        msg = f'{name+":":<20} {t_tot:.2e}'
        if global_tot is not None:
            msg += f' ({100*t_tot/global_tot:.1f} %)'
        if count > 1:
            msg += f' [{count}, {t_min:.2e} / {t_avg:.2e} / {t_max:.2e}]'
        print('\t'*tabs+f'{msg}')
    else:
        if debug:
            print('!'*60)
            print(f'Timing for "{tag}" not found')
            print('!'*60)

    return t_tot

def readProfile(name, max_lvl):
    # Get file handle
    f = h5py.File(name, 'r')
    print(name)

    info = f['info']
    print('git commit: ', info['git-commit'][()][0])

    ranks = info['ranks'][()][0]
    print('ranks: ', ranks)

    db = f['timings']

    print('Timing format: min / avg / max')
    print('#'*40)

    # Walltime
    indent = 0
    printTiming('Walltime', 'Walltime', db, indent)
    print('\n')

    # Init
    indent = 0
    printHeader('Initialization', indent)
    printTiming('Model', 'createSimulation', db, indent+1)
    printTiming('Simulation', 'Simulation::preRun', db, indent+1)

    # Computation
    indent = 0
    tot = printTiming('Computation', 'Simulation::mainRun', db, indent)
    printTiming('evolve', 'Pseudospectral::Coordinator::evolve', db, indent+1, global_tot = tot)
    printTiming('explicitEquations', 'Pseudospectral::Coordinator::explicitEquations', db, indent+2, global_tot = tot)
    if max_lvl > 1:
        printTiming('trivial', 'Pseudospectral::Coordinator::explicitEquations-trivial', db, indent+3, global_tot = tot)
        printTiming('diagnostic', 'Pseudospectral::Coordinator::explicitEquations-diagnostic', db, indent+3, global_tot = tot)
        printTiming('prognostic', 'Pseudospectral::Coordinator::explicitEquations-prognostic', db, indent+3, global_tot = tot)
    printTiming('computeNonlinear', 'Pseudospectral::Coordinator::computeNonlinear', db, indent+2, global_tot = tot)
    printTiming('updatePhysical', 'Pseudospectral::Coordinator::updatePhysical', db, indent+3, global_tot = tot)
    if max_lvl > 1:
        printTiming('prepareSpectral', 'Transform::BackwardConfigurator::prepareSpectral', db, indent+4, global_tot = tot)
        printTiming('project1D', 'Transform::BackwardConfigurator::project1D', db, indent+4, global_tot = tot)
        if max_lvl > 2:
            printTiming('pre', 'Transform::BackwardConfigurator::project1D-pre', db, indent+5, global_tot = tot)
            printTiming('transform', 'Transform::BackwardConfigurator::project1D-transform', db, indent+5, global_tot = tot)
            printTiming('post', 'Transform::BackwardConfigurator::project1D-post', db, indent+5, global_tot = tot)
        printTiming('project2D', 'Transform::BackwardConfigurator::project2D', db, indent+4, global_tot = tot)
        if max_lvl > 2:
            printTiming('pre', 'Transform::BackwardConfigurator::project2D-pre', db, indent+5, global_tot = tot)
            printTiming('transform', 'Transform::BackwardConfigurator::project2D-transform', db, indent+5, global_tot = tot)
            printTiming('post', 'Transform::BackwardConfigurator::project2D-post', db, indent+5, global_tot = tot)
        printTiming('preparePhysical', 'Transform::BackwardConfigurator::preparePhysical', db, indent+4, global_tot = tot)
        printTiming('projectND', 'Transform::BackwardConfigurator::projectND', db, indent+4, global_tot = tot)
        if max_lvl > 2:
            printTiming('pre', 'Transform::BackwardConfigurator::projectND-pre', db, indent+5, global_tot = tot)
            printTiming('transform', 'Transform::BackwardConfigurator::projectND-transform', db, indent+5, global_tot = tot)
            printTiming('post', 'Transform::BackwardConfigurator::projectND-post', db, indent+5, global_tot = tot)
    tag = 'Transform::ForwardConfigurator::nonlinearTerm'
    printTiming('nonlinearTerm', tag, db, indent+3, global_tot = tot)
    for t in db:
        if t.startswith(tag + '-'):
            printTiming(t.split(tag + '-')[1], t, db, indent+4, global_tot = tot)
    printTiming('updateSpectral', 'Pseudospectral::Coordinator::updateSpectral', db, indent+3, global_tot = tot)
    if max_lvl > 1:
        printTiming('integrateND', 'Transform::ForwardConfigurator::integrateND', db, indent+4, global_tot = tot)
        if max_lvl > 2:
            printTiming('pre', 'Transform::ForwardConfigurator::integrateND-pre', db, indent+5, global_tot = tot)
            printTiming('transform', 'Transform::ForwardConfigurator::integrateND-transform', db, indent+5, global_tot = tot)
            printTiming('post', 'Transform::ForwardConfigurator::integrateND-post', db, indent+5, global_tot = tot)
        printTiming('integrate2D', 'Transform::ForwardConfigurator::integrate2D', db, indent+4, global_tot = tot)
        if max_lvl > 2:
            printTiming('pre', 'Transform::ForwardConfigurator::integrate2D-pre', db, indent+5, global_tot = tot)
            printTiming('transform', 'Transform::ForwardConfigurator::integrate2D-transform', db, indent+5, global_tot = tot)
            printTiming('post', 'Transform::ForwardConfigurator::integrate2D-post', db, indent+5, global_tot = tot)
        printTiming('integrate1D', 'Transform::ForwardConfigurator::integrate1D', db, indent+4, global_tot = tot)
        if max_lvl > 2:
            printTiming('pre', 'Transform::ForwardConfigurator::integrate1D-pre', db, indent+5, global_tot = tot)
            printTiming('transform', 'Transform::ForwardConfigurator::integrate1D-transform', db, indent+5, global_tot = tot)
            printTiming('post', 'Transform::ForwardConfigurator::integrate1D-post', db, indent+5, global_tot = tot)
        printTiming('updateEquation', 'Transform::ForwardConfigurator::updateEquation', db, indent+4, global_tot = tot)
    printTiming('solveEquations', 'Pseudospectral::Coordinator::solveEquations', db, indent+2, global_tot = tot)
    if max_lvl > 1:
        printTiming('trivial-before', 'Pseudospectral::Coordinator::solveEquations-trivialBefore', db, indent+3, global_tot = tot)
        printTiming('diagnostic-before', 'Pseudospectral::Coordinator::solveEquations-diagnosticBefore', db, indent+3, global_tot = tot)
        printTiming('prognostic', 'Pseudospectral::Coordinator::solveEquations-prognostic', db, indent+3, global_tot = tot)
        printTiming('diagnostic-after', 'Pseudospectral::Coordinator::solveEquations-diagnosticAfter', db, indent+3, global_tot = tot)
        printTiming('trivial-after', 'Pseudospectral::Coordinator::solveEquations-trivialAfter', db, indent+3, global_tot = tot)
    tag = 'Pseudospectral::Coordinator::updateEquations'
    printTiming('updateEquations', tag, db, indent+2, global_tot = tot)
    for t in db:
        if t.startswith(tag + '-'):
            printTiming(t.split(tag + '-')[1], t, db, indent+4, global_tot = tot)
    printTiming('finalizeTimestep', 'Pseudospectral::Coordinator::finalizeTimestep', db, indent+2, global_tot = tot)
    printTiming('IO', 'Simulation::writeOutput', db, indent+1, global_tot = tot)
    printTiming('Stats', 'SimulationIoControl::writeStats', db, indent+2, global_tot = tot)
    printTiming('updateHeavyAscii', 'SimulationIoTools::updateHeavyAscii', db, indent+2, global_tot = tot)
    printTiming('Ascii', 'SimulationIoControl::writeAscii', db, indent+2, global_tot = tot)
    printTiming('Hdf5', 'SimulationIoControl::writeHdf5', db, indent+2, global_tot = tot)
    printTiming('Diagnostics', 'Pseudospectral::Coordinator::writeDiagnostics', db, indent+2, global_tot = tot)

    # Cleanup
    indent = 0
    printTiming('Cleanup', 'Simulation::postRun', db, indent)

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        # Path
        path = '.'
        lvl = 3
    elif(len(sys.argv) < 3):
        # Path
        path = sys.argv[1]
        lvl = 3
    else:
        path = sys.argv[1]
        lvl = int(sys.argv[2])

    name = '/profile.hdf5'
    readProfile(path+name, lvl)
