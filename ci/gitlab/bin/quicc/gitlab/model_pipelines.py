#
# Base and derived classes to generate yml model pipelines for QuICC
#

from quicc.gitlab.models import default_configs
from quicc.gitlab.pipelines import config
from quicc.gitlab.yaml import base_yaml

"""Base class, defines the steps to build the models"""
class base_model(base_yaml):
    def __init__(self, cnfs):
        self.lib_path = '${CI_CACHE_FOLDER}/quicc-${QUICC_VERSION_TAG}'
        # this contains all model configs
        self.cnfs = cnfs
        self.config = {}
        # pipeline actions
        self.actions = [ self.model_build_yaml ]

    def pipe(self):
        for c in self.cnfs:
            self.config_model(c)
            for a in self.actions:
                a()

    def set_file_name(self):
        self.file_name = '.quicc_models'

    def config_model(self, cnf):
        self.model_name = cnf.name
        self.tag = cnf.tag
        self.test = cnf.test
        self.config_name = cnf.shortname()
        self.pipe_line_name = cnf.fullname()
        self.full_tag = cnf.fulltag()
        self.variant_extension = cnf.variantextension()
        self.model_path = self.lib_path+'-'+self.config_name

    def model_build_yaml(self):
        self.config['.'+self.pipe_line_name] = {
            'stage': 'model-build-and-test',
            'script':
            [
                '/QuICC.src/ci/gitlab/bin/mpi_lock.sh hostname',
                # first install the library to the shared folder
                # note, only one model per pipeline does this step
                'cd /QuICC.src/build',
                '/QuICC.src/ci/gitlab/bin/mpi_lock.sh \"cmake ..' \
                    ' -DCMAKE_INSTALL_PREFIX='+self.lib_path+ \
                    ' -DQUICC_TESTSUITE_POLYNOMIAL=OFF' \
                    ' -DQUICC_TESTSUITE_FRAMEWORK=OFF' \
                    ' -DQUICC_TESTSUITE_TRANSFORM=OFF' \
                    ' -DQUICC_TESTSUITE_PROFILING=OFF' \
                    ' -DQUICC_TESTSUITE_SPARSESM=OFF' \
                    ' && make install\"',
                # then move to the shared location to build the model
                'export QUICC_ROOT='+self.model_path,
                'echo $QUICC_ROOT',
                'mkdir -p $QUICC_ROOT/build',
                'cd $QUICC_ROOT/build',
                # echo command is added to make lock unique per each pipeline
                '/QuICC.src/ci/gitlab/bin/mpi_lock.sh \"cmake /QuICC.src --log-level=VERBOSE -DQUICC_MODEL='+self.model_name+' -DQUICC_TESTSUITE_MODEL=ON -DQUICC_GITHUB_PROTOCOL=https -DQUICC_USE_SYSTEM_QUICC=ON -Dquicc_DIR='+self.lib_path+'/share/quicc/cmake && bash -c \'time make -j $(grep processor /proc/cpuinfo | wc -l) && echo ${QUICC_VERSION_TAG}_'+self.config_name+'\'\"',
            ],
        }

"""Extend base class to test the models"""
class test_model(base_model):
    def __init__(self, cnf):
        super(test_model, self).__init__(cnf)
        self.actions.extend([self.model_test_yaml])

    def model_test_yaml(self):
        # check before actually adding the test
        if self.test:
            self.config['.'+self.pipe_line_name]['script'].extend(
                [
                    # export manually non-default path
                    'export PYTHONPATH=$QUICC_ROOT/build/lib/python:$PYTHONPATH',
                    # cannot run the mpi model through ctest within a srun command
                    'cd '+self.model_path+'/build/Models/'+self.model_name+'/TestSuite/Benchmarks/_data/'+self.full_tag,
                    self.model_path+'/build/Models/'+self.model_name+'/Executables/'+self.config_name+'Model',
                    'cd '+self.model_path+'/build',
                    '/QuICC.src/ci/gitlab/bin/mpi_lock.sh \"ctest --no-tests=error --verbose -R ValidateBenchmark'+self.config_name+'Model'+self.variant_extension+'$ && echo ${QUICC_VERSION_TAG}\"'
                ]
            )

"""Extend base class to time the models"""
class perf_model(base_model):
    def __init__(self, cnf):
        super(perf_model, self).__init__(cnf)
        self.actions.extend([self.model_time_yaml])

    def set_file_name(self):
        self.file_name = '.quicc_models_perf'

    def model_time_yaml(self):
        self.config['.'+self.pipe_line_name]['script'].extend(
            [
                # export manually non-default path
                'export PYTHONPATH=$QUICC_ROOT/build/lib/python:$PYTHONPATH',
                # use reframe to setup stage folder only
                '/QuICC.src/ci/gitlab/bin/mpi_lock.sh reframe ' \
                ' -c /QuICC.src/ci/reframe/quicc.py -r -v ' \
                ' -S quicc_root=$QUICC_ROOT ' \
                ' -S target_executable=ls ' \
                ' -S model='+self.model_name+ \
                ' -S tag='+self.tag+ \
                ' --keep-stage-files --skip-sanity-check --skip-performance-check --force-local',
                # run reframe only on one rank
                '/QuICC.src/ci/reframe/bin/mpi_wrapper.sh ' \
                ' "reframe -c /QuICC.src/ci/reframe/quicc.py -r -v'\
                ' -S quicc_root=$QUICC_ROOT'\
                ' -S model='+self.model_name+ \
                ' -S tag='+self.tag+ \
                ' --performance-report --keep-stage-files --force-local --dont-restage"' \
                ' "cd stage/generic/default/builtin/run_quicc_63_127_127_1_1_5; ' \
                +self.model_path+'/build/Models/'+self.model_name+'/Executables/'+self.config_name+'Model"'
            ]
        )
