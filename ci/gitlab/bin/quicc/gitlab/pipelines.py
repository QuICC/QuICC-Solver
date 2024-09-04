#
# Base and derived classes to generate yml pipelines for QuICC
#
import hashlib
from typing import NamedTuple
from quicc.gitlab.models import default_configs
from quicc.gitlab.yaml import base_yaml

"""Config struct to generate a pipeline"""
class config(NamedTuple):
    tag: str
    backend: str

backend2nodeSize = {
    "daint-mc": 72,
    "daint-gpu": 24,
    "alps-a100": 64
}

"""Base class, defines a pipeline that build the docker image, the library and cleans up the runner"""
class base_pipeline(base_yaml):
    def __init__(self, cnf):
        image_location = '$CSCS_REGISTRY_PATH'
        self.cs_base_yml = 'https://gitlab.com/cscs-ci/recipes/-/raw/master/templates/v2/.ci-ext.yml'
        self.backend = cnf.backend
        self.base_docker = 'ci/docker/baseimage/Dockerfile_quicc_baseimage_'+self.backend
        self.cpus_full_node = backend2nodeSize[self.backend]
        base_md5sum = hashlib.md5(open(self.base_docker[3:], 'rb').read()).hexdigest()
        self.tag = cnf.tag
        image = 'quicc_'+cnf.tag+'_'+cnf.backend+':$CI_COMMIT_SHA'
        self.path_image = image_location+'/'+image
        self.base_path_image = f'{image_location}/baseimage/quicc_baseimage_{cnf.backend}:{base_md5sum}'
        self.docker = 'ci/docker/Dockerfile_'+cnf.tag

        # pipeline actions
        self.actions = [ self.base_yaml ]

    """append cleanup to extend more easily the pipelines"""
    def pipe(self):
        for a in self.actions:
            a()
        self.cleanup_yaml()

    def set_file_name(self):
        self.file_name = '.quicc_'+self.tag+'_'+self.backend

    def base_yaml(self):
        self.config = {
            'include':
                [
                    {'remote': self.cs_base_yml},
                    '/ci/gitlab/.cscs_builders.yml',
                    '/ci/gitlab/.cscs_runners.yml',
                ],
            'stages':
                [
                    'build_base',
                    'build', # build stages are running on the Kubernetes cluster
                ],
            'build-quicc-base':
                {
                    'extends': '.'+self.backend+'_builder',
                    'stage': 'build_base',
                    'variables':
                        {
                            'DOCKERFILE': self.base_docker,
                            'PERSIST_IMAGE_NAME': self.base_path_image,
                        },
                },
            'build-quicc':
                {
                    'extends': '.'+self.backend+'_builder',
                    'stage': 'build',
                    'variables':
                        {
                            'DOCKERFILE': self.docker,
                            'PERSIST_IMAGE_NAME': self.path_image,
                            'DOCKER_BUILD_ARGS': f"""["BASEIMAGE={self.base_path_image}"]"""
                        },
                },
            }

    def cleanup_yaml(self):
        self.config['stages'].extend(
                [
                    'cleanup',
                ],
            )
        if (self.backend == 'gpu'):
            runner = '.container-runner-daint-gpu'
        else:
            runner = '.container-runner-daint'
        self.config['ci-cache-cleanup'] = {
            'extends': runner,
            'stage': 'cleanup',
            'image': self.path_image,
            'script':
                [
                    'rm -Rf $CI_CACHE_FOLDER/*'
                ],
            }

"""Add library testing to the base pipeline"""
class libtest_pipeline(base_pipeline):
    def __init__(self, cnf):
        super(libtest_pipeline, self).__init__(cnf)
        self.actions.extend([self.testlib_yaml])

    def testlib_yaml(self):
        self.config['include'].extend(
                [
                    '/ci/gitlab/.quicc_tests.yml',
                ],
            )
        self.config['stages'].extend(
                [
                    'test', # test stage is running on PizDaint (on 1 node)
                ],
            )
        self.config['test-quicc-lib'] = {
                'extends':
                    [
                        '.test-lib',
                        '.'+self.backend+'_runner'
                    ],
                'image': self.path_image,
            }
        if (self.tag == 'mpi'):
            tasks = '4'
            self.config['test-quicc-mpi-lib'] = {
                'extends':
                    [
                        '.test-lib-mpi',
                        '.'+self.backend+'_runner'
                    ],
                'image': self.path_image,
                'variables':
                    {
                        'SLURM_NTASKS': tasks,
                        'SLURM_NTASKS_PER_NODE': tasks,
                        'SLURM_CPUS_PER_TASK': str(self.cpus_full_node//int(tasks))
                    },
            }

"""Add library timing to the libtest pipeline"""
class libtime_pipeline(libtest_pipeline):
    def __init__(self, cnf):
        super(libtime_pipeline, self).__init__(cnf)
        self.actions.extend([self.timelib_yaml])

    def timelib_yaml(self):
        self.config['time-quicc-lib'] = {
                'extends':
                    [
                        '.time-lib-'+self.backend,
                        '.'+self.backend+'_runner'
                    ],
                'image': self.path_image,
            }

"""Add library sweep timing to the base pipeline and changes default name of the yml file"""
class libtime_sweep_pipeline(base_pipeline):
    def __init__(self, cnf):
        super(libtime_sweep_pipeline, self).__init__(cnf)
        self.actions.extend([self.timelib_yaml])

    def set_file_name(self):
        self.file_name = '.quicc_'+self.tag+'_'+self.backend+'_perf'

    def timelib_yaml(self):
        self.config['include'].extend(
                [
                    '/ci/gitlab/.quicc_tests.yml',
                ],
            )
        self.config['stages'].extend(
                [
                    'test',
                ],
            )
        self.config['time-quicc-lib'] = {
                'extends':
                    [
                        '.time-lib-sweep-'+self.backend,
                        '.'+self.backend+'_runner'
                    ],
                'image': self.path_image,
            }


"""Add model testing to the libtime pipeline"""
class model_pipeline(libtime_pipeline):
    def __init__(self, cnf):
        super(model_pipeline, self).__init__(cnf)
        self.actions.extend([self.model_yaml])

    def model_yaml(self):
        self.config['include'].extend(
                [
                    '/ci/gitlab/.quicc_models.yml',
                ],
            )
        self.config['stages'].extend(
                [
                    'model-build-and-test',
                ],
            )
        for mode_config in default_configs(self.tag):
            model = mode_config.fullname()
            tasks = str(mode_config.tasks)
            self.config[model] = {
                    'extends':
                        [
                            '.'+model,
                            '.'+self.backend+'_runner'
                        ],
                    'image': self.path_image,
                    'variables':
                    {
                        # the image is pulled in the lib test stage
                        'PULL_IMAGE': 'NO',
                        'SLURM_NTASKS': tasks,
                        'SLURM_NTASKS_PER_NODE': tasks,
                        'SLURM_CPUS_PER_TASK': str(self.cpus_full_node//int(tasks)),
                        'QUICC_VERSION_TAG': self.tag
                    },
                }

"""Add model testing to the libtest pipeline"""
class model_pipeline_notiming(libtest_pipeline):
    def __init__(self, cnf):
        super(model_pipeline_notiming, self).__init__(cnf)
        self.actions.extend([self.model_yaml])

    def model_yaml(self):
        self.config['include'].extend(
                [
                    '/ci/gitlab/.quicc_models.yml',
                ],
            )
        self.config['stages'].extend(
                [
                    'model-build-and-test',
                ],
            )
        for mode_config in default_configs(self.tag):
            model = mode_config.fullname()
            tasks = str(mode_config.tasks)
            self.config[model] = {
                    'extends':
                        [
                            '.'+model,
                            '.'+self.backend+'_runner'
                        ],
                    'image': self.path_image,
                    'variables':
                    {
                        # the image is pulled in the lib test stage
                        'PULL_IMAGE': 'NO',
                        'SLURM_NTASKS': tasks,
                        'SLURM_NTASKS_PER_NODE': tasks,
                        'SLURM_CPUS_PER_TASK': str(self.cpus_full_node//int(tasks)),
                        'QUICC_VERSION_TAG': self.tag
                    },
                }


"""Add model timing to the base pipeline and changes default name of the yml file"""
class model_perf_pipeline(base_pipeline):
    def __init__(self, cnf):
        super(model_perf_pipeline, self).__init__(cnf)
        self.actions.extend([self.model_yaml])

    def set_file_name(self):
        self.file_name = '.quicc_'+self.tag+'_'+self.backend+'_models_perf'

    def model_yaml(self):
        self.config['include'].extend(
                [
                    '/ci/gitlab/.quicc_models_perf.yml',
                ],
            )
        self.config['stages'].extend(
                [
                    'model-build-and-test',
                ],
            )
        for mode_config in default_configs('perf'):
            model = mode_config.fullname()
            if mode_config.tasks == -1:
                tasks = str(int(self.cpus_full_node/2))
            else:
                tasks = str(mode_config.tasks)
            self.config[model] = {
                    'extends':
                        [
                            '.'+model,
                            '.'+self.backend+'_runner'
                        ],
                    'image': self.path_image,
                    'variables':
                    {
                        'SLURM_NTASKS': tasks,
                        'SLURM_NTASKS_PER_NODE': tasks,
                        'SLURM_CPUS_PER_TASK': str(self.cpus_full_node//int(tasks)),
                        'QUICC_VERSION_TAG': self.tag
                    },
                }
