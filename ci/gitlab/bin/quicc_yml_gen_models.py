#
# Script to generate yml files for QuICC model tests
#
from quicc.gitlab.models import default_configs
from quicc.gitlab.model_pipelines import perf_model, test_model

if __name__ == '__main__':
    # Test pipeline
    # print(default_configs('serial')+default_configs('mpi'))
    pipe = test_model(default_configs('serial')+default_configs('mpi'))
    pipe.write()

    # Perf pipeline
    pipe = perf_model(default_configs('perf'))
    pipe.write()
