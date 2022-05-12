#
# Script to generate yml files for QuICC
#
import os
import sys
import yaml
from typing import NamedTuple

class config(NamedTuple):
    tag: str
    backend: str
    model: bool

def populateYaml(cnf):
    imageLocation = 'art.cscs.ch/contbuild/testing/anfink'
    image = 'quicc_'+cnf.tag+':latest'
    configIn = {
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
                        'DOCKERFILE': 'docker/Dockerfile_'+cnf.tag,
                        'PERSIST_IMAGE_NAME': image
                    },
            },
        'include':
            [
                '/ci/.daint_runner.yml',
                '/ci/.quicc_tests.yml',
                '/ci/.quicc_models.yml'
            ],
        'test-quicc-lib_'+cnf.tag:
            {
                'extends':
                    [
                        '.test-lib',
                        '.'+cnf.backend
                    ],
                'image': imageLocation+'/'+image
            },

        }
    if(cnf.model):
        if cnf.tag == 'mpi':
            tasks = '4'
        else:
            tasks = '1'
        for model in ['BoussinesqSphereDynamo', 'BoussinesqShellDynamo',
            'BoussinesqSphereRTC', 'BoussinesqShellRTC',
            'BoussinesqSphereTC', 'BoussinesqShellTC']:
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
        yml = validateYaml(populateYaml(c))
        gitlab_yml = open('.quicc_'+c.tag+'_'+c.backend+'.yml','w')
        gitlab_yml.write('# This file is generated, do not modify!\n')
        gitlab_yml.write(yml)
        gitlab_yml.close()
