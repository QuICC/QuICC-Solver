"""@file quicc.py"""

# Copyright 2016-2022 Swiss National Supercomputing Centre (CSCS/ETH Zurich)

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.core.runtime as rt
import contextlib

@rfm.simple_test
class run_quicc(rfm.RunOnlyRegressionTest):
    """QuICC model performance test
    """
    dims = parameter(['63;127;127'])
    quicc_root = variable(str, value='./')
    exe_path = ''
    use_tool = variable(str, value='notool')  # scorep+p,perftools,extrae,mpip
    cfg = variable(str, value='parameters.cfg')
    rpt = variable(str, value='OUT_stdout')
    model = variable(str, value='notset')
    tag = variable(str, value='notset')
    target_executable = variable(str, value='notset')
    csv_rpt = variable(str, value='rpt.csv')
    compute_nodes = parameter([1])
    gpus_per_cn = parameter([1])
    steps = parameter([5])
    cscs_systems = ['daint:mc', 'manali:mc', 'dom:gpu', 'dom:mc']
    local_systems = ['generic', 'g-X1']
    valid_systems = cscs_systems + local_systems
    valid_prog_environs = ['PrgEnv-cray', 'PrgEnv-gnu', 'builtin']
    # only PrgEnv-cray supports affinity correctly
    use_multithreading = False
    system = rt.runtime().system

    refs =  {   'BoussinesqSphereDynamoExplicitModel': {
                    1: {
                        'icelake': {
                            'rpt_t_Initialisation': (34, -0.1, 1.0, 's'),
                            'rpt_t_PreRun': (4.4, -0.1, 1.0, 's'),
                            'rpt_t_Computation': (74.0, -0.5, 1.0, 's'),
                            'rpt_t_PostRun': (0.05, -1.0, 2.0, 's'),
                            'rpt_t_Total': (120, -1.0, 1.0, 's'),
                        },
                        'broadwell': {
                            'rpt_t_Initialisation': (12.8, -0.05, 0.4, 's'),
                            'rpt_t_PreRun': (2.2, -0.10, 1.0, 's'),
                            'rpt_t_Computation': (10.0, -0.05, 0.4, 's'),
                            'rpt_t_PostRun': (0.2, -0.5, 1.2, 's'),
                        }
                    }
                }
            }

    @run_after('init')
    def set_paths(self):
        """Set various paths """
        self.sourcesdir = self.quicc_root+'/build/Models/'+self.model+'/TestSuite/Benchmarks/_refdata/'+self.tag
        self.exe_path = self.quicc_root+'/build/Models/'+self.model+'/Executables'

    @run_after('init')
    def set_target(self):
        """Set target """
        # `ls` exe is needed to trick reframe into setting up the stage directory
        # without running the model
        if(self.target_executable != 'ls'):
            self.target_executable = self.model+self.tag+'Model'

    @run_before('run')
    def set_num_task(self):
        """Set num tasks based on machine"""
        proc = self.current_partition.processor
        self.num_tasks_per_node = proc.num_cores
        self.num_tasks = self.compute_nodes * self.num_tasks_per_node
        self.time_limit = '20m'

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
        if self.use_tool == 'perftools':
            self.executable = 'pat_run'
            self.executable_name = f'{self.exe_path}/{self.target_executable}'
            self.executable_opts = ['-r', self.executable_name]
        else:
            self.executable = f'{self.exe_path}/{self.target_executable}'

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

    @run_before('run')
    def set_cfg(self):
        nr = int(self.dims.split(";")[0])
        nl = int(self.dims.split(";")[1])
        nm = int(self.dims.split(";")[2])
        self.prerun_cmds += [
            '# --- sed ---',
            rf'sed -i "s@<cpus>[0-9]\+</cpus>@<cpus>{self.num_tasks}</cpus>@" '
            f'{self.cfg}',
            rf'sed -i "s@<dim1D>[0-9]\+</dim1D>@<dim1D>{nr}</dim1D>@" '
            f'{self.cfg}',
            rf'sed -i "s@<dim2D>[0-9]\+</dim2D>@<dim2D>{nl}</dim2D>@" '
            f'{self.cfg}',
            rf'sed -i "s@<dim3D>[0-9]\+</dim3D>@<dim3D>{nm}</dim3D>@" '
            f'{self.cfg}',
            rf'sed -i "s@<sim>-[0-9]\+</sim>@<sim>-{self.steps}</sim>@" '
            f'{self.cfg}',
            #
            rf'sed -i "s@algorithm>[a-z0-9]\+</algorithm>@algorithm>tubular</algorithm>@" {self.cfg}',
            rf'sed -i "s@<ascii>[0-9]\+</ascii>@<ascii>0</ascii>@" {self.cfg}',
            rf'sed -i "s@<hdf5>[0-9]\+</hdf5>@<hdf5>0</hdf5>@" {self.cfg}',
            rf'sed -i "s@<enable>1</enable>@<enable>0</enable>@" {self.cfg}',
        ]

    @sanity_function
    def assert_sanity(self):
        regex1 = r'Total \[s\]'
        sanity_list = [
            sn.assert_found(regex1, self.rpt),
        ]
        if self.use_tool == 'scorep+p':
            sanity_list.append(
                sn.assert_found(r'Estimated aggregate size of event trace',
                                self.rpt_score)
            )
        elif self.use_tool == 'extrae':
            sanity_list.append(
                sn.assert_found(r'Congratulations! \S+ has been generated.',
                                self.stdout)
            )
        elif self.use_tool == 'mpip':
            sanity_list.append(
                sn.assert_found(r'^mpiP: Storing mpiP output in',
                                self.stdout)
            )

        return sn.all(sanity_list)

    def report_params(self, regex_str):
        ''' check timers '''
        if regex_str not in ('dim1D', 'dim2D', 'dim3D', 'cpus', 'Timesteps'):
            raise ValueError(f'illegal value in argument ({regex_str!r})')

        # rpt = 'OUT_stdout'
        regex = r'^\s+%s: (\d+)' % regex_str
        return sn.extractsingle(regex, self.rpt, 1, int)

    @performance_function('')
    def rpt_dim1(self):
        return self.report_params('dim1D')

    @performance_function('')
    def rpt_dim2(self):
        return self.report_params('dim2D')

    @performance_function('')
    def rpt_dim3(self):
        return self.report_params('dim3D')

    @performance_function('')
    def rpt_cpus(self):
        return self.report_params('cpus')

    @performance_function('')
    def rpt_steps(self):
        return self.report_params('Timesteps')


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
                pname: self.refs[self.target_executable][self.compute_nodes][arch]
            }

    def report_time(self, regex_str):
        """
        *********** Execution time information ***********
        --------------------------------------------------
           Initialisation [s]: 14.161 / 14.169 / 14.155
                PreRun [s]: 3.627 / 3.631 / 3.624
              Computation [s]: 2.575 / 2.575 / 2.575
                PostRun [s]: 3.427 / 3.427 / 3.427
        """
        if regex_str not in ('Initialisation', 'PreRun', 'Computation',
                             'PostRun', 'Total'):
            raise ValueError(f'illegal value in argument ({regex_str!r})')

        # rpt = 'OUT_stdout'
        regex = r'^\s+%s \[s\]: (\S+)' % regex_str
        return sn.extractsingle(regex, self.rpt, 1, float)

    @performance_function('s')
    def rpt_t_Initialisation(self):
        return self.report_time('Initialisation')

    @performance_function('s')
    def rpt_t_PreRun(self):
        return self.report_time('PreRun')

    @performance_function('s')
    def rpt_t_Computation(self):
        return self.report_time('Computation')

    @performance_function('s')
    def rpt_t_PostRun(self):
        return self.report_time('PostRun')

    @performance_function('s')
    def rpt_t_Total(self):
        return self.report_time('Total')
