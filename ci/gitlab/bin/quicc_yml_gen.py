#
# Script to generate yml files for QuICC
#
import os
import sys
import yaml
from typing import NamedTuple
from quicc_defaults import defaultConfigs

class config(NamedTuple):
    tag: str
    backend: str
    model: bool

def populateYaml(cnf):
    imageLocation = '$CSCS_REGISTRY_PATH'
    baseImage = 'ubuntu-quicc-buildbase:22.04'
    image = 'quicc_'+cnf.tag+':latest'
    csBaseYml = 'https://gitlab.com/cscs-ci/recipes/-/raw/master/templates/v2/.cscs.yml'
    configIn = {
        'include':
            [
                {'remote': csBaseYml},
                '/ci/gitlab/.daint_runner.yml',
                '/ci/gitlab/.quicc_tests.yml',
                '/ci/gitlab/.quicc_models.yml',
            ],
        'stages':
            [
                'build', # build stage is running on the Kubernetes cluster
                'test', # test stage is running on PizDaint (on 1 node)
                'model-build-and-test',
                'cleanup',
            ],
        'build-quicc_'+cnf.tag:
            {
                'tags':
                    [
                        'docker_jfrog'
                    ],
                'stage': 'build',
                'script':
                    [
                        'true'
                    ],
                'variables':
                    {
                        'DOCKERFILE': 'ci/docker/Dockerfile_'+cnf.tag,
                        'PERSIST_IMAGE_NAME': imageLocation+'/'+image,
                        'DOCKER_BUILD_ARGS': f"""["BASEIMAGE={imageLocation}/{baseImage}"]"""
                    },
            },
        'test-quicc-lib_'+cnf.tag:
            {
                'extends':
                    [
                        '.test-lib',
                        '.'+cnf.backend
                    ],
                'image': imageLocation+'/'+image
            },
        'time-quicc-lib_'+cnf.tag:
            {
                'extends':
                    [
                        '.time-lib',
                        '.'+cnf.backend
                    ],
                'image': imageLocation+'/'+image
            },
        'ci-cache-cleanup':
            {
                'extends': '.daint',
                'stage': 'cleanup',
                'image': imageLocation+'/'+image,
                'script':
                    [
                        'rm -Rf $CI_CACHE_FOLDER/*'
                    ],
            },
        }
    if(cnf.model):
        for modCnf in defaultConfigs():
            if cnf.tag in modCnf.ignore:
                continue
            model = modCnf.fullname()
            if cnf.tag == 'mpi':
                tasks = str(modCnf.tasks)
            else:
                tasks = '1'
            configIn[model+'_'+cnf.tag] = {
                    'extends':
                        [
                            '.'+model,
                            '.'+cnf.backend
                        ],
                    'image': imageLocation+'/'+image,
                    'variables':
                    {
                        # the image is pulled in the lib test stage
                        'PULL_IMAGE': 'NO',
                        'SLURM_NTASKS': tasks,
                        'SLURM_NTASKS_PER_NODE': tasks,
                        'SLURM_CPUS_PER_TASK': str(72//int(tasks)),
                        'QUICC_VERSION_TAG': cnf.tag
                    },
                }

    configOut = yaml.dump(configIn, default_flow_style=False,
        sort_keys=False)
    return configOut

def validateYaml(yml):
    try:
        yaml.safe_load(yml)
        return yml
    except:
        sys.exit('Failed to validate yml.')

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        confs = [
            config('serial', 'cpu', True),
            config('mp', 'cpu', False),
            config('mpi', 'cpu', True),
            ]
    else:
        confs = config(sys.argv[1], sys.argv[2], sys.argv[3])

    for c in confs:
        # yml = validateYaml(populateYaml(c))
        yml = populateYaml(c)
        gitlab_yml = open('.quicc_'+c.tag+'_'+c.backend+'.yml','w')
        gitlab_yml.write('# This file is generated, do not modify!\n')
        gitlab_yml.write(yml)
        gitlab_yml.close()
