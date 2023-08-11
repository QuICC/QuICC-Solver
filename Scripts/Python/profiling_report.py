#! /usr/bin/env python3

import h5py
import subprocess
import sys
import numpy as np

class Formatter:
    def __init__(self, name, skip = 0, debug = False):
        # Get file handle
        f = h5py.File(name, 'r')
        print(name)

        self.info = f['info']
        commit = self.info['git-commit'][()][0]
        print('git commit: ', commit)

        self.ranks = self.info['ranks'][()][0]
        print('ranks: ', self.ranks)

        self.db = f['timings']

        self.skip = skip

        self.debug = debug

    def printHeader(self, name, tabs):
        print('\t'*tabs+f'{name+":":<20}')

    def printTiming(self, name, tag, tabs, global_tot = None):
        t_tot = None
        if tag in self.db:
            t0 = 0
            timings = self.db[tag]
            count = timings['count'][()]
            sze = timings['time'][()].shape[0]
            ts = np.reshape(timings['time'][()],(self.ranks,sze//self.ranks)).T
            if count > 1:
                t0 = self.skip
            ts = ts[t0:,:]
            count -= t0
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
            if self.debug:
                print('!'*60)
                print(f'Timing for "{tag}" not found')
                print('!'*60)

        return t_tot

def readProfile(name, max_lvl):
    f = Formatter(name, skip=0, debug = True)

    print('Timing format: min / avg / max')
    print('#'*40)

    # Walltime
    indent = 0
    f.printTiming('Walltime', 'Walltime', indent)
    print('\n')

    # Init
    indent = 0
    f.printHeader('Initialization', indent)
    f.printTiming('Model', 'createSimulation', indent+1)
    f.printTiming('Simulation', 'Simulation::preRun', indent+1)

    # Computation
    indent = 0
    tot = f.printTiming('Computation', 'Simulation::mainRun', indent)
    f.printTiming('evolve', 'Pseudospectral::Coordinator::evolve', indent+1, global_tot = tot)
    f.printTiming('explicitEquations', 'Pseudospectral::Coordinator::explicitEquations', indent+2, global_tot = tot)
    if max_lvl > 1:
        f.printTiming('trivial', 'Pseudospectral::Coordinator::explicitEquations-trivial', indent+3, global_tot = tot)
        f.printTiming('diagnostic', 'Pseudospectral::Coordinator::explicitEquations-diagnostic', indent+3, global_tot = tot)
        f.printTiming('prognostic', 'Pseudospectral::Coordinator::explicitEquations-prognostic', indent+3, global_tot = tot)
    f.printTiming('computeNonlinear', 'Pseudospectral::Coordinator::computeNonlinear', indent+2, global_tot = tot)
    f.printTiming('updatePhysical', 'Pseudospectral::Coordinator::updatePhysical', indent+3, global_tot = tot)
    if max_lvl > 1:
        f.printTiming('prepareSpectral', 'Transform::BackwardConfigurator::prepareSpectral', indent+4, global_tot = tot)
        f.printTiming('project1D', 'Transform::BackwardConfigurator::project1D', indent+4, global_tot = tot)
        if max_lvl > 2:
            f.printTiming('pre', 'Transform::BackwardConfigurator::project1D-pre', indent+5, global_tot = tot)
            f.printTiming('transform', 'Transform::BackwardConfigurator::project1D-transform', indent+5, global_tot = tot)
            f.printTiming('post', 'Transform::BackwardConfigurator::project1D-post', indent+5, global_tot = tot)
        f.printTiming('project2D', 'Transform::BackwardConfigurator::project2D', indent+4, global_tot = tot)
        if max_lvl > 2:
            f.printTiming('pre', 'Transform::BackwardConfigurator::project2D-pre', indent+5, global_tot = tot)
            f.printTiming('transform', 'Transform::BackwardConfigurator::project2D-transform', indent+5, global_tot = tot)
            f.printTiming('post', 'Transform::BackwardConfigurator::project2D-post', indent+5, global_tot = tot)
        f.printTiming('preparePhysical', 'Transform::BackwardConfigurator::preparePhysical', indent+4, global_tot = tot)
        f.printTiming('projectND', 'Transform::BackwardConfigurator::projectND', indent+4, global_tot = tot)
        if max_lvl > 2:
            f.printTiming('pre', 'Transform::BackwardConfigurator::projectND-pre', indent+5, global_tot = tot)
            f.printTiming('transform', 'Transform::BackwardConfigurator::projectND-transform', indent+5, global_tot = tot)
            f.printTiming('post', 'Transform::BackwardConfigurator::projectND-post', indent+5, global_tot = tot)
    tag = 'Transform::ForwardConfigurator::nonlinearTerm'
    f.printTiming('nonlinearTerm', tag, indent+3, global_tot = tot)
    for t in f.db:
        if t.startswith(tag + '-'):
            f.printTiming(t.split(tag + '-')[1], t, indent+4, global_tot = tot)
    f.printTiming('updateSpectral', 'Pseudospectral::Coordinator::updateSpectral', indent+3, global_tot = tot)
    if max_lvl > 1:
        f.printTiming('integrateND', 'Transform::ForwardConfigurator::integrateND', indent+4, global_tot = tot)
        if max_lvl > 2:
            f.printTiming('pre', 'Transform::ForwardConfigurator::integrateND-pre', indent+5, global_tot = tot)
            f.printTiming('transform', 'Transform::ForwardConfigurator::integrateND-transform', indent+5, global_tot = tot)
            f.printTiming('post', 'Transform::ForwardConfigurator::integrateND-post', indent+5, global_tot = tot)
        f.printTiming('integrate2D', 'Transform::ForwardConfigurator::integrate2D', indent+4, global_tot = tot)
        if max_lvl > 2:
            f.printTiming('pre', 'Transform::ForwardConfigurator::integrate2D-pre', indent+5, global_tot = tot)
            f.printTiming('transform', 'Transform::ForwardConfigurator::integrate2D-transform', indent+5, global_tot = tot)
            f.printTiming('post', 'Transform::ForwardConfigurator::integrate2D-post', indent+5, global_tot = tot)
        f.printTiming('integrate1D', 'Transform::ForwardConfigurator::integrate1D', indent+4, global_tot = tot)
        if max_lvl > 2:
            f.printTiming('pre', 'Transform::ForwardConfigurator::integrate1D-pre', indent+5, global_tot = tot)
            f.printTiming('transform', 'Transform::ForwardConfigurator::integrate1D-transform', indent+5, global_tot = tot)
            f.printTiming('post', 'Transform::ForwardConfigurator::integrate1D-post', indent+5, global_tot = tot)
        f.printTiming('updateEquation', 'Transform::ForwardConfigurator::updateEquation', indent+4, global_tot = tot)
    f.printTiming('solveEquations', 'Pseudospectral::Coordinator::solveEquations', indent+2, global_tot = tot)
    if max_lvl > 1:
        f.printTiming('trivial-before', 'Pseudospectral::Coordinator::solveEquations-trivialBefore', indent+3, global_tot = tot)
        f.printTiming('diagnostic-before', 'Pseudospectral::Coordinator::solveEquations-diagnosticBefore', indent+3, global_tot = tot)
        f.printTiming('prognostic', 'Pseudospectral::Coordinator::solveEquations-prognostic', indent+3, global_tot = tot)
        f.printTiming('diagnostic-after', 'Pseudospectral::Coordinator::solveEquations-diagnosticAfter', indent+3, global_tot = tot)
        f.printTiming('trivial-after', 'Pseudospectral::Coordinator::solveEquations-trivialAfter', indent+3, global_tot = tot)
    tag = 'Pseudospectral::Coordinator::updateEquations'
    f.printTiming('updateEquations', tag, indent+2, global_tot = tot)
    for t in f.db:
        if t.startswith(tag + '-'):
            f.printTiming(t.split(tag + '-')[1], t, indent+4, global_tot = tot)
    f.printTiming('finalizeTimestep', 'Pseudospectral::Coordinator::finalizeTimestep', indent+2, global_tot = tot)
    f.printTiming('IO', 'Simulation::writeOutput', indent+1, global_tot = tot)
    f.printTiming('Stats', 'SimulationIoControl::writeStats', indent+2, global_tot = tot)
    f.printTiming('updateHeavyAscii', 'SimulationIoTools::updateHeavyAscii', indent+2, global_tot = tot)
    f.printTiming('Ascii', 'SimulationIoControl::writeAscii', indent+2, global_tot = tot)
    f.printTiming('Hdf5', 'SimulationIoControl::writeHdf5', indent+2, global_tot = tot)
    f.printTiming('Diagnostics', 'Pseudospectral::Coordinator::writeDiagnostics', indent+2, global_tot = tot)

    # Cleanup
    indent = 0
    f.printTiming('Cleanup', 'Simulation::postRun', indent)

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
