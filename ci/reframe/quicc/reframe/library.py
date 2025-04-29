"""@file quicc_library.py"""

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.core.runtime as rt
import contextlib
import quicc.reframe.utils as quicc_utils
from collections import defaultdict

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

class testBase(rfm.RunOnlyRegressionTest):
    """QuICC transform performance test
    """
    quicc_root = variable(str, value='./')
    exe_path = ''
    target_executable = variable(str, value='ctest')
    test = ''
    region = []
    csv_rpt = variable(str, value='rpt.csv')
    cscs_systems = ['daint:mc', 'manali:mc', 'dom:gpu', 'dom:mc', 'alps:a100', 'alps:gh200']
    local_systems = ['generic', 'g-X1', 'ph-fangorn']
    valid_systems = cscs_systems + local_systems
    valid_prog_environs = ['PrgEnv-cray', 'PrgEnv-gnu', 'builtin']
    is_serial_test = True
    # only PrgEnv-cray supports affinity correctly
    use_multithreading = False
    system = rt.runtime().system
    refs =  {}

    @run_before('run')
    def set_num_task(self):
        """Set num tasks based on machine"""
        proc = self.current_partition.processor
        self.num_tasks_per_node = proc.num_cores
        self.num_tasks = 1
        self.time_limit = '30m'

    @run_before('run')
    def set_prgenv(self):
        if(self.system.name in self.cscs_systems):
            self.modules = [
                'cray-fftw', 'cray-hdf5-parallel', 'cray-python', 'cray-tpsl',
                'Boost'
            ]
        else:
            self.modules = []

    @run_before('run')
    def set_exe(self):
        args = '--test-dir '+self.quicc_root+'/build -V -R'
        self.executable = f'{self.target_executable} {args} {self.test}'

        self.omp = {
            'dom': {'omp': '12', '-c': '24'},
            'manali': {'omp': '16', '-c': '32'},
        }
        if(self.system.name in self.cscs_systems):
            self.job.launcher.options = [
                '-n', str(self.num_tasks), '--cpu-bind=cores'
            ]

            self.prerun_cmds += [
                'module list',
                'echo SLURM_NTASKS=$SLURM_NTASKS',
                'echo SLURM_NTASKS_PER_NODE=$SLURM_NTASKS_PER_NODE',
                'echo SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK',
                'echo SLURM_DISTRIBUTION=$SLURM_DISTRIBUTION',
                'echo starttime=`date +%s`',
            ]
            self.postrun_cmds += [
                'echo stoptime=`date +%s`',
                'echo "job done"',
                'echo "SLURMD_NODENAME=$SLURMD_NODENAME"',
                'echo "SLURM_JOBID=$SLURM_JOBID"',
                # f'# dims={self.dims} ',
                # 'rm -f core*',
            ]
        elif(not self.is_serial_test):
            self.job.launcher.options = [
                '-n', str(self.num_tasks_per_node)
            ]

    @run_before('run')
    def set_perf_reference(self):
        proc = self.current_partition.processor
        device = self.current_partition.devices
        pname = self.current_partition.fullname
        arch = proc.arch

        # fix for gpu systems
        if pname in ('daint:gpu', 'dom:gpu'):
            arch = 'p100'
        elif pname in ('alps:a100'):
            arch = 'a100'
        elif pname in ('alps:gh200'):
            arch = 'gh200'
        # fix for CI
        elif pname == 'generic:default':
            if(arch == 'haswell'):
                arch = 'p100'
            elif(arch == 'zen3'):
                arch = 'a100'
            elif(arch == 'neoverse_v2'):
                arch = 'gh200'

        with contextlib.suppress(KeyError):
            self.reference = {
                pname: self.refs[arch]
            }

    @sanity_function
    def assert_sanity(self):
        regex1 = r'All tests passed'
        sanity_list = [
            sn.assert_found(regex1, self.stdout),
        ]

        return sn.all(sanity_list)

    def report_time(self, region, group):
        """
47: IALegendreProjector::applyOperators 100 times
47:      stats by rank                         min             max             avg
47:      high watermark memory (MB):       29.3929     29.3929     29.3929
47:      memory delta (MB)         :             0           0           0
47:      time (s)                  :      0.001117     0.00185  0.00130549
        """
        regex = r"{}".format(region)+r'[\s\S]*?time \(s\) [\D]*(\S+)[\D]*(\S+)[\D]*(\S+)'
        # third group is the average value
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
            for id, l1 in db[test].items():
                for split, l2 in l1.items():
                    for arch, l3 in l2.items():
                        # Group by base region to collect all values
                        for region, time in l3.items():
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

                        # Add references with adjusted times using the full region name
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

class testStateFile(testBase):
    """QuICC transform performance test
    """

    @performance_function('s')
    def perfIoAvg(self):
        return self.report_time_avg(self.region)

    @performance_function('s')
    def perfIoMin(self):
        return self.report_time_min(self.region)

    @performance_function('s')
    def perfIoMax(self):
        return self.report_time_max(self.region)

    @performance_function('s')
    def perfIoScalarsAvg(self):
        return self.report_time_avg(self.region+'-scalars')

    @performance_function('s')
    def perfIoScalarsMin(self):
        return self.report_time_min(self.region+'-scalars')

    @performance_function('s')
    def perfIoScalarsMax(self):
        return self.report_time_max(self.region+'-scalars')

    @performance_function('s')
    def perfIoVectorsAvg(self):
        return self.report_time_avg(self.region+'-vectors')

    @performance_function('s')
    def perfIoVectorsMin(self):
        return self.report_time_min(self.region+'-vectors')

    @performance_function('s')
    def perfIoVectorsMax(self):
        return self.report_time_max(self.region+'-vectors')
