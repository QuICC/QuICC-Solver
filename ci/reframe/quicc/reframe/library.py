"""@file quicc_library.py"""

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.core.runtime as rt
import contextlib

class testBase(rfm.RunOnlyRegressionTest):
    """QuICC transform performance test
    """
    quicc_root = variable(str, value='./')
    exe_path = ''
    target_executable = variable(str, value='ctest')
    test = ''
    region = ''
    csv_rpt = variable(str, value='rpt.csv')
    steps = 100
    cscs_systems = ['daint:mc', 'manali:mc', 'dom:gpu', 'dom:mc']
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
        pname = self.current_partition.fullname
        arch = proc.arch

        # fix for gpu systems
        if pname in ('daint:gpu', 'dom:gpu'):
            arch = 'p100'
        # fix for CI
        elif pname == 'generic:default':
            if(arch == 'haswell'):
                arch = 'p100'

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

        regexSteps = str(self.steps)+r' times'
        sanity_list.append(sn.assert_found(regexSteps, self.stdout))

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

        # third group is average value
        return sn.extractsingle(regex, self.stdout, group, float)

    def report_time_avg(self, region):
        return self.report_time(region, 3)

    def report_time_min(self, region):
        return self.report_time(region, 1)

    def report_time_max(self, region):
        return self.report_time(region, 2)

class testTransform(testBase):
    """QuICC transform performance test
    """

    @performance_function('s')
    def applyOperatorsAvg(self):
        return self.report_time_avg(self.region)

    @performance_function('s')
    def applyOperatorsMin(self):
        return self.report_time_min(self.region)

    @performance_function('s')
    def applyOperatorsMax(self):
        return self.report_time_max(self.region)

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
