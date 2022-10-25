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
    csv_rpt = variable(str, value='rpt.csv')
    steps = 100
    cscs_systems = ['daint:mc', 'manali:mc', 'dom:gpu', 'dom:mc']
    local_systems = ['generic', 'g-X1']
    valid_systems = cscs_systems + local_systems
    valid_prog_environs = ['PrgEnv-cray', 'PrgEnv-gnu', 'builtin']
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
        if(self.system.name in ['daint', 'manali', 'dom']):
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

    @run_before('run')
    def set_perf_reference(self):
        proc = self.current_partition.processor
        pname = self.current_partition.fullname
        if pname in ('daint:gpu', 'dom:gpu'):
            arch = 'p100'
        else:
            arch = proc.arch

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
        regex = r'applyOperators[\s\S]*time \(s\) [\D]*(\S+)[\D]*(\S+)[\D]*(\S+)'

        # third group is average value
        return sn.extractsingle(regex, self.stdout, group, float)

    def report_time_avg(self, region):
        return self.report_time(region, 3)

    def report_time_min(self, region):
        return self.report_time(region, 1)

    def report_time_max(self, region):
        return self.report_time(region, 2)


    @performance_function('s')
    def applyOperatorsAvg(self):
        return self.report_time_avg('applyOperators')

    @performance_function('s')
    def applyOperatorsMin(self):
        return self.report_time_min('applyOperators')

    @performance_function('s')
    def applyOperatorsMax(self):
        return self.report_time_max('applyOperators')


# actual tests
@rfm.simple_test
class testALegendreTests_Poly_P_projector_id206(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id206_ulp19000'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00012, -0.15, 0.3, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000138504, -0.15, 0.2, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id207(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id207_ulp180000'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.0013, -0.15, 0.3, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000935358, -0.15, 0.2, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split96_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id108_ulp20800_split96_0'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00930405, -0.1, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0100948, -0.1, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split288_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id108_ulp20800_split288_0'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00367698, -0.1, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00395184, -0.1, 0.1, 's'),
                },
            }




