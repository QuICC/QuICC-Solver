#
# Script to generate yml files for QuICC pipelines on CSCS CI
#
import os
from typing import NamedTuple
from quicc.gitlab.pipelines import config, libtest_pipeline, \
    model_pipeline, model_pipeline_notiming, model_perf_pipeline, \
    libtime_sweep_pipeline

if __name__ == '__main__':
    # Base pipelines
    base_confs = [
        config('mp', 'daint-mc')
        ]
    for c in base_confs:
        pipe = libtest_pipeline(c)
        pipe.write()

    # Model without timing pipelines
    model_notiming_confs = [
        config('serial', 'daint-mc')
        ]
    for c in model_notiming_confs:
        pipe = model_pipeline_notiming(c)
        pipe.write()

    # Model and Timing pipelines
    model_confs = [
        config('mpi', 'daint-mc'),
        config('kk', 'daint-mc'),
        config('kkgpu', 'daint-gpu'),
        config('kkgpu', 'alps-a100')
        ]
    for c in model_confs:
        pipe = model_pipeline(c)
        pipe.write()

    # Library perf pipelines
    lib_perf_confs = [
        config('mpi', 'daint-mc'),
        config('kkgpu', 'daint-gpu'),
        config('kkgpu', 'alps-a100')
    ]
    for c in lib_perf_confs:
        pipe = libtime_sweep_pipeline(c)
        pipe.write()

    # Model perf pipelines
    model_perf_confs = [
        config('mpi', 'daint-mc'),
    ]
    for c in model_perf_confs:
        pipe = model_perf_pipeline(c)
        pipe.write()
