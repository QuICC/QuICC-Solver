#
# Script to generate yml files for QuICC model tests
#
import os
import sys
import yaml
from quicc_defaults import config, defaultConfigs

def populateYaml(cnf):
    libPath = '${CI_CACHE_FOLDER}/quicc-${QUICC_VERSION_TAG}'
    modelPath = libPath+'-'+cnf.name+cnf.tag
    configIn = { '.'+cnf.name+cnf.tag:
        {
            'stage': 'model-build-and-test',
            'script':
            [
                'hostname',
                # first install the library to the shared folder
                # note, only one model per pipeline does this step
                'cd /QuICC.src/build',
                '/QuICC.src/ci/gitlab/bin/mpi_lock.sh \"cmake ..' \
                    ' -DCMAKE_INSTALL_PREFIX='+libPath+ \
                    ' -DQUICC_TESTSUITE_POLYNOMIAL=OFF' \
                    ' -DQUICC_TESTSUITE_FRAMEWORK=OFF' \
                    ' -DQUICC_TESTSUITE_TRANSFORM=OFF' \
                    ' -DQUICC_TESTSUITE_PROFILING=OFF' \
                    ' -DQUICC_TESTSUITE_SPARSESM=OFF' \
                    ' && make install\"',
                # then move to the shared location to build the model
                'mkdir -p '+modelPath+'/build',
                'cd '+modelPath+'/build',
                # echo command is added to make lock unique per each pipeline
                '/QuICC.src/ci/gitlab/bin/mpi_lock.sh \"cmake /QuICC.src --log-level=VERBOSE -DQUICC_MODEL='+cnf.name+' -DQUICC_TESTSUITE_MODEL=ON -DQUICC_GITHUB_PROTOCOL=https -DQUICC_USE_SYSTEM_QUICC=ON -Dquicc_DIR='+libPath+'/lib/cmake/quicc && bash -c \'time make -j $(grep processor /proc/cpuinfo | wc -l) && echo ${QUICC_VERSION_TAG}_'+cnf.name+cnf.tag+'\'\"',

            ],
        }
    }
    if cnf.test:
        configIn['.'+cnf.name+cnf.tag]['script'].extend(
            [
                # export manually non-default path
                'export PYTHONPATH='+modelPath+'/build/lib/python',
                # cannot run the mpi model through ctest within a srun command
                'cd '+modelPath+'/build/Models/'+cnf.name+'/TestSuite/Benchmarks/_data/'+cnf.tag,
                modelPath+'/build/Models/'+cnf.name+'/Executables/'+cnf.name+cnf.tag+'Model',
                'cd '+modelPath+'/build',
                '/QuICC.src/ci/gitlab/bin/mpi_lock.sh \"ctest --no-tests=error --verbose -R ValidateBenchmark'+cnf.name+cnf.tag+' && echo ${QUICC_VERSION_TAG}\"'
            ]
        )

    configOut = yaml.dump(configIn, default_flow_style=False,
        sort_keys=False, width=200)
    return configOut

def validateYaml(yml):
    try:
        yaml.safe_load(yml)
        return yml
    except:
        sys.exit('Failed to validate yml.')

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        confs = defaultConfigs()
    else:
        confs = config(sys.argv[1], sys.argv[2], sys.argv[3])

    gitlab_yml = open('.quicc_models.yml','w')
    gitlab_yml.write('# This file is generated, do not modify!\n')
    for c in confs:
        yml = validateYaml(populateYaml(c))
        gitlab_yml.write(yml)
    gitlab_yml.close()
