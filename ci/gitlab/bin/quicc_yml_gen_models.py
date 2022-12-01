#
# Script to generate yml files for QuICC model tests
#
from quicc.gitlab.models import default_configs
from quicc.gitlab.model_pipelines import perf_model, test_model

if __name__ == '__main__':
    # Test pipeline
    model_test_configs = []
    for t in ['serial', 'mpi', 'kk', 'kkgpu']:
        model_test_configs.extend(default_configs(t))
    pipe = test_model(model_test_configs)
    pipe.write()

    # Perf pipeline
    pipe = perf_model(default_configs('perf'))
    pipe.write()
