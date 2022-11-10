#
# Script to generate yml files for QuICC
#
import os
from typing import NamedTuple
from quicc.gitlab.pipelines import config, libtest_pipeline, \
    model_pipeline, perf_pipeline

if __name__ == '__main__':
    # Base pipelines
    base_confs = [
        config('mp', 'cpu'),
        ]
    for c in base_confs:
        pipe = libtest_pipeline(c)
        pipe.write()

    # Model pipelines
    model_confs = [
        config('serial', 'cpu'),
        config('mpi', 'cpu'),
        ]
    for c in model_confs:
        pipe = model_pipeline(c)
        pipe.write()

    # Perf pipelines
    perf_confs = [
        config('mpi', 'cpu'),
    ]
    for c in perf_confs:
        pipe = perf_pipeline(c)
        pipe.write()
