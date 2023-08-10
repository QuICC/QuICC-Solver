#! /usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import subprocess
import sys
import numpy as np

def printHeader(name, tabs):
    print('\t'*tabs+f'{name+":":<20}')

def printTiming(name, tag, db, tabs, debug = True):
    if tag in db:
        timings = db[tag]
        count = timings['count'][()]
        ts = timings['time'][()]
        t_tot = np.average(ts)*count
        t_min = np.min(ts)
        t_avg = np.average(ts)
        t_max = np.max(ts)
        if count > 1:
            print('\t'*tabs+f'{name+":":<20} {t_tot:.2e} ({count}, {t_min:.2e} / {t_avg:.2e} / {t_max:.2e})')
        else:
            print('\t'*tabs+f'{name+":":<20} {t_tot:.2e}')
    else:
        if debug:
            print('!'*60)
            print(f'Timing for "{tag}" not found')
            print('!'*60)

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
    printTiming('Computation', 'Simulation::mainRun', db, indent)
    printTiming('evolve', 'Pseudospectral::Coordinator::evolve', db, indent+1)
    printTiming('explicitEquations', 'Pseudospectral::Coordinator::explicitEquations', db, indent+2)
    if max_lvl > 1:
        printTiming('trivial', 'Pseudospectral::Coordinator::explicitEquations-trivial', db, indent+3)
        printTiming('diagnostic', 'Pseudospectral::Coordinator::explicitEquations-diagnostic', db, indent+3)
        printTiming('prognostic', 'Pseudospectral::Coordinator::explicitEquations-prognostic', db, indent+3)
    printTiming('computeNonlinear', 'Pseudospectral::Coordinator::computeNonlinear', db, indent+2)
    printTiming('updatePhysical', 'Pseudospectral::Coordinator::updatePhysical', db, indent+3)
    if max_lvl > 1:
        printTiming('prepareSpectral', 'Transform::BackwardConfigurator::prepareSpectral', db, indent+4)
        printTiming('project1D', 'Transform::BackwardConfigurator::project1D', db, indent+4)
        if max_lvl > 2:
            printTiming('pre', 'Transform::BackwardConfigurator::project1D-pre', db, indent+5)
            printTiming('transform', 'Transform::BackwardConfigurator::project1D-transform', db, indent+5)
            printTiming('post', 'Transform::BackwardConfigurator::project1D-post', db, indent+5)
        printTiming('project2D', 'Transform::BackwardConfigurator::project2D', db, indent+4)
        if max_lvl > 2:
            printTiming('pre', 'Transform::BackwardConfigurator::project2D-pre', db, indent+5)
            printTiming('transform', 'Transform::BackwardConfigurator::project2D-transform', db, indent+5)
            printTiming('post', 'Transform::BackwardConfigurator::project2D-post', db, indent+5)
        printTiming('preparePhysical', 'Transform::BackwardConfigurator::preparePhysical', db, indent+4)
        printTiming('projectND', 'Transform::BackwardConfigurator::projectND', db, indent+4)
        if max_lvl > 2:
            printTiming('pre', 'Transform::BackwardConfigurator::projectND-pre', db, indent+5)
            printTiming('transform', 'Transform::BackwardConfigurator::projectND-transform', db, indent+5)
            printTiming('post', 'Transform::BackwardConfigurator::projectND-post', db, indent+5)
    printTiming('nonlinearTerm', 'Transform::ForwardConfigurator::nonlinearTerm', db, indent+3)
    printTiming('updateSpectral', 'Pseudospectral::Coordinator::updateSpectral', db, indent+3)
    if max_lvl > 1:
        printTiming('integrateND', 'Transform::ForwardConfigurator::integrateND', db, indent+4)
        if max_lvl > 2:
            printTiming('pre', 'Transform::ForwardConfigurator::integrateND-pre', db, indent+5)
            printTiming('transform', 'Transform::ForwardConfigurator::integrateND-transform', db, indent+5)
            printTiming('post', 'Transform::ForwardConfigurator::integrateND-post', db, indent+5)
        printTiming('integrate2D', 'Transform::ForwardConfigurator::integrate2D', db, indent+4)
        if max_lvl > 2:
            printTiming('pre', 'Transform::ForwardConfigurator::integrate2D-pre', db, indent+5)
            printTiming('transform', 'Transform::ForwardConfigurator::integrate2D-transform', db, indent+5)
            printTiming('post', 'Transform::ForwardConfigurator::integrate2D-post', db, indent+5)
        printTiming('integrate1D', 'Transform::ForwardConfigurator::integrate1D', db, indent+4)
        if max_lvl > 2:
            printTiming('pre', 'Transform::ForwardConfigurator::integrate1D-pre', db, indent+5)
            printTiming('transform', 'Transform::ForwardConfigurator::integrate1D-transform', db, indent+5)
            printTiming('post', 'Transform::ForwardConfigurator::integrate1D-post', db, indent+5)
        printTiming('updateEquation', 'Transform::ForwardConfigurator::updateEquation', db, indent+4)
    printTiming('solveEquations', 'Pseudospectral::Coordinator::solveEquations', db, indent+2)
    if max_lvl > 1:
        printTiming('trivial-before', 'Pseudospectral::Coordinator::solveEquations-trivialBefore', db, indent+3)
        printTiming('diagnostic-before', 'Pseudospectral::Coordinator::solveEquations-diagnosticBefore', db, indent+3)
        printTiming('prognostic', 'Pseudospectral::Coordinator::solveEquations-prognostic', db, indent+3)
        printTiming('diagnostic-after', 'Pseudospectral::Coordinator::solveEquations-diagnosticAfter', db, indent+3)
        printTiming('trivial-after', 'Pseudospectral::Coordinator::solveEquations-trivialAfter', db, indent+3)
    printTiming('updateEquations', 'Pseudospectral::Coordinator::updateEquations', db, indent+2)
    printTiming('finalizeTimestep', 'Pseudospectral::Coordinator::finalizeTimestep', db, indent+2)
    printTiming('IO', 'Simulation::writeOutput', db, indent+1)
    printTiming('Stats', 'SimulationIoControl::writeStats', db, indent+2)
    printTiming('updateHeavyAscii', 'SimulationIoTools::updateHeavyAscii', db, indent+2)
    printTiming('Ascii', 'SimulationIoControl::writeAscii', db, indent+2)
    printTiming('Hdf5', 'SimulationIoControl::writeHdf5', db, indent+2)
    printTiming('Diagnostics', 'Pseudospectral::Coordinator::writeDiagnostics', db, indent+2)

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
