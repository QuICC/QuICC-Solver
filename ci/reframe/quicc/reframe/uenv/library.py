# Copyright 2024 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os
import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.udeps as udeps
import quicc.reframe.utils as quicc_utils
from collections import defaultdict
import uenv

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

class quicc_download(rfm.RunOnlyRegressionTest):
    descr = 'Fetch QuICC sources code'
    sourcesdir = None
    executable = 'git'
    executable_opts = [
        'clone', 'git@github.com:QuICC/QuICC.git'
        '&&', 'cd', 'QuICC'
        '&&', 'git', 'checkout', 'reframe_uenv'
    ]
    local = True

    @sanity_function
    def validate_download(self):
        return sn.assert_eq(self.job.exitcode, 0)

class quicc_build(rfm.CompileOnlyRegressionTest):
    descr = 'Build QuICC'
    valid_systems = ['+uenv']
    valid_prog_environs = ['+mpi']
    build_system = 'CMake'
    sourcedir = None
    maintainers = ['dganellari']
    quicc_sources = fixture(quicc_download, scope='session')
    # NOTE: required so that the build stage is performed on
    # a compute node using an sbatch job.
    # This will force the uenv and view to be loaded using
    # "#SBATCH --uenv=" etc
    build_locally = False

    @run_before('compile')
    def prepare_build(self):
        self.uarch = uenv.uarch(self.current_partition)
        self.build_system.builddir = os.path.join(self.stagedir, 'build')

        self.prebuild_cmds = [
            f'rsync -a {self.quicc_sources.stagedir}/QuICC/ .'
        ]

        self.build_system.config_opts = [
            ' -DCMAKE_BUILD_TYPE=Release',
            ' -DQUICC_USE_KOKKOS=ON',
            ' -DQUICC_EIGEN_ENABLE_VECTORIZATION=ON',
            ' -DQUICC_PROFILE_LEVEL=3',
            ' -DQUICC_PROFILE_BACKEND=native',
            ' -DQUICC_PROFILE_NATIVE_SAMPLE=500',
            ' -DQUICC_TESTSUITE_TYPES=ON',
            ' -DQUICC_TESTSUITE_TRANSFORM=ON',
            ' -DQUICC_GITHUB_PROTOCOL=https',
        ]
        # set architecture-specific flags
        if self.uarch == 'gh200':
            self.build_system.config_opts += [
                ' -DCMAKE_CUDA_COMPILER=nvcc',
                ' -DCMAKE_CUDA_ARCHITECTURES=90',
                ' -DQUICC_USE_KOKKOS_CUDA=ON',
                ' -DQUICC_EIGEN_ENABLE_CUDA=ON',
                ' -DQUICC_USE_VKFFT=ON'
                ' -DQUICC_USE_CUFFT=ON',
                ' -DCMAKE_CUDA_FLAGS_RELEASE="-expt-extended-lambda"',
            ]
        elif self.uarch == 'a100':
            self.build_system.config_opts += [
                ' -DCMAKE_CUDA_COMPILER=nvcc',
                ' -DCMAKE_CUDA_ARCHITECTURES=80',
                ' -DQUICC_USE_KOKKOS_CUDA=ON',
                ' -DQUICC_EIGEN_ENABLE_CUDA=ON',
                ' -DQUICC_USE_VKFFT=ON',
                ' -DQUICC_USE_CUFFT=ON',
                ' -DCMAKE_CUDA_FLAGS_RELEASE="-expt-extended-lambda"',
            ]
        elif self.uarch == 'zen2':
            self.build_system.config_opts += []

        self.build_system.max_concurrency = 64

        self.build_system.make_opts = []

class testBase(rfm.RunOnlyRegressionTest):
    valid_systems = ['+uenv']
    valid_prog_environs = ['+mpi']
    target_executable = variable(str, value='ctest')
    maintainers = ['dganellari']

    test = ''
    region = []
    is_serial_test = True
    refs = {}

    quicc_build = fixture(quicc_build, scope='environment')

    @run_before('run')
    def set_num_tasks(self):
        """Set num tasks based on machine"""
        proc = self.current_partition.processor
        self.num_tasks_per_node = proc.num_cores
        self.num_tasks = 1
        self.time_limit = '30m'

    @run_before('run')
    def set_exec(self):
        self.executable = self.target_executable
        self.executable_opts = [
            '--test-dir', f'{self.quicc_build.stagedir}/build', '-V', '-R', self.test
        ]

    @run_before('run')
    def set_perf_reference(self):
        self.uarch = uenv.uarch(self.current_partition)

        if (self.uarch is not None) and (self.uarch in self.refs):
            self.reference = {
                self.current_partition.fullname:
                    self.refs[self.uarch]
            }

    @sanity_function
    def assert_sanity(self):
        return sn.assert_found(r'All tests passed', self.stdout)

    def report_time(self, region, group):
        """
47: IALegendreProjector::applyOperators 100 times
47:      stats by rank                         min             max             avg
47:      high watermark memory (MB):       29.3929     29.3929     29.3929
47:      memory delta (MB)         :             0           0           0
47:      time (s)                  :      0.001117     0.00185  0.00130549
        """
        regex = r"{}".format(region)+r'[\s\S]*?time \(s\) [\D]*(\S+)[\D]*(\S+)[\D]*(\S+)'

        # third group is average value
        return sn.extractsingle(regex, self.stdout, group, float)

    def report_iterations(self, region):
        """Extract number of iterations
47: IALegendreProjector::applyOperators 100 times
        """
        regex = r"{}[\s\S]*?(\d+) times".format(region)
        return sn.extractsingle(regex, self.stdout, 1, int)

    def report_time_avg(self, region):
        return self.report_time(region, 3)

    def report_time_min(self, region):
        return self.report_time(region, 1)

    def report_time_max(self, region):
        return self.report_time(region, 2)

    def report_iter(self, region):
        return self.report_iterations(region)

    def report_adjusted_time(self, region):
        """Compute and report adjusted timing using (avg * n_iter - max)/(n_iter - 1)"""
        avg = self.report_time_avg(region)
        max_time = self.report_time_max(region)
        iterations = self.report_iter(region)
        return self.compute_adjusted_time(avg, max_time, iterations)

    @run_before('performance')
    def set_perf_variables(self):
        """Build the dictionary with all the performance variables."""

        for r in self.region:
            self.perf_variables[r+'Avg'] = sn.make_performance_function(self.report_time_avg(r), 's')
            self.perf_variables[r+'Min'] = sn.make_performance_function(self.report_time_min(r), 's')
            self.perf_variables[r+'Max'] = sn.make_performance_function(self.report_time_max(r), 's')
            self.perf_variables[r+'Iter'] = sn.make_performance_function(self.report_iter(r), '')
            self.perf_variables[r+'AdjAvg'] = sn.make_performance_function(
                lambda r=r: self.report_adjusted_time(r), 's'
            )

class testTransform(testBase):
    """QuICC transform performance test
    """

class splitTestTransform(testTransform):
    """QuICC transform performance split test
    """

    splitting = variable(int, value=8*12)
    testId = variable(int, value=108)

    def __init__(self):
        self.refs_db = nested_dict(3,dict)
        self.init_references()

    def init_references(self):
        pass

    def add_reference(self, tId, split, arch, region, timings):
        """ Add new timing reference
        """

        self.refs_db[tId][split][arch][region] = timings

    def compute_adjusted_time(self, avg: float, max_time: float, iterations: int) -> float:
        """Compute adjusted timing using (avg * n_iter - max)/(n_iter - 1)"""
        if iterations <= 1:
            return avg
        return (avg * iterations - max_time)/(iterations - 1)

    def read_references(self, db_file, test):
        """Read references and compute adjusted time from JSON data"""
        db = quicc_utils.read_reference_timings(db_file)
        temp_data = {}

        # Extract test data from DB if available
        if test in db:
            for id,l1 in db[test].items():
                for split,l2 in l1.items():
                    for arch,l3 in l2.items():
                        # Group by base region to collect all values
                        for region,time in l3.items():
                            base_region = region.rsplit('Avg')[0].rsplit('Max')[0].rsplit('Iter')[0]

                            if base_region not in temp_data or temp_data[base_region]['arch'] != arch:
                                temp_data[base_region] = {
                                    'avg': None, 'max': None, 'iter': None,
                                    'id': id, 'split': split, 'arch': arch
                                }

                            if region.endswith('Avg'):
                                temp_data[base_region]['avg'] = float(time)
                            elif region.endswith('Max'):
                                temp_data[base_region]['max'] = float(time)
                            elif region.endswith('Iter'):
                                temp_data[base_region]['iter'] = int(time)

                        # Determine tolerance values based on test ID
                        if 'kokkos' in test:
                            lower_bound = -0.4
                            upper_bound = 0.3
                        else:
                            lower_bound = -0.25
                            upper_bound = 0.2

                        # Add references with adjusted times
                        for base_region, values in temp_data.items():
                            if all(values.get(k) is not None for k in ['avg', 'max', 'iter']):
                                adjusted_time = self.compute_adjusted_time(
                                    values['avg'], values['max'], values['iter']
                                )
                                self.add_reference(
                                    int(values['id']),
                                    int(values['split']),
                                    values['arch'],
                                    base_region + 'AdjAvg',  # Use the full region name
                                    (float(adjusted_time), lower_bound, upper_bound, 's')
                                )

    @run_before('compile')
    def setupTest(self):
        self.test = f'prof_{self.__class__.__name__}_id{self.testId}_ulp.*_split{self.splitting}_0'
        self.refs = self.refs_db[self.testId][self.splitting]

