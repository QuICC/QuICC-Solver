# Daint
## Model
```sh
export QUICC_ROOT=</path/to/QuICC>
reframe -c $QUICC_ROOT/ci/reframe/quicc.py --keep-stage-files --failure-stats --performance-report --system daint:mc -r -p PrgEnv-gnu -S use_tool=none -S quicc_root=$QUICC_ROOT
```
## Test
```sh
export QUICC_ROOT=</path/to/QuICC>
reframe -c $QUICC_ROOT/ci/reframe/quicc_library.py -r -S quicc_root=$QUICC_ROOT --keep-stage-files --performance-report --force-local
```

# Local machine
## Model
```sh
export QUICC_ROOT=</path/to/QuICC>
reframe -C $QUICC_ROOT/ci/reframe/settings.py -c $QUICC_ROOT/ci/reframe/quicc.py -r -v -S quicc_root=$QUICC_ROOT --performance-report
```
## Test
```sh
export QUICC_ROOT=</path/to/QuICC>
export PYTHONPATH=$QUICC_ROOT/build/lib/python:$PYTHONPATH
export PYTHONPATH=$QUICC_ROOT/ci/reframe:$PYTHONPATH
reframe -C $QUICC_ROOT/ci/reframe/settings.py -c $QUICC_ROOT/ci/reframe/quicc_library_[cpu|gpu].py -r -S quicc_root=$QUICC_ROOT --performance-report --force-local --exec-policy=serial
```

# CI
## Model
```sh
export QUICC_ROOT=</path/to/QuICC>
reframe -c $QUICC_ROOT/ci/reframe/quicc.py -r -S quicc_root=$QUICC_ROOT
 -S target_executable=BoussinesqSphereDynamoExplicitModel
--performance-report --force-local
```
## Test
```sh
export QUICC_ROOT=</path/to/QuICC>
reframe -c $QUICC_ROOT/ci/reframe/quicc_library_[cpu|gpu].py -r -S quicc_root=$QUICC_ROOT --performance-report --force-local --exec-policy=serial
```

# Extracting timings from performance report and update references
- Download report from CI as text file, say quicc_time_cpu.txt for arch=broadwell
```sh
export QUICC_ROOT=</path/to/QuICC>
export PYTHONPATH=$QUICC_ROOT/ci/reframe:$PYTHONPATH
python -c "import quicc.reframe.utils as utils;d = utils.extract_timings('quicc_time_cpu.txt', 'broadwell');utils.write_timings('new_daint_mc_cpu.json', d)"
python -c "import quicc.reframe.utils as utils;d = utils.update_timings('new_daint_mc_cpu.json', 'cpu.json')"
```

# Plot default data
- Download report from CI as text file or use existing json files
```sh
export QUICC_ROOT=</path/to/QuICC>
export PYTHONPATH=$QUICC_ROOT/ci/reframe:$PYTHONPATH
python -c "import quicc.reframe.utils as utils;utils.make_default_plots(file_cpu = 'quicc_time_cpu.txt', file_gpu = 'quicc_time_gpu.txt', save = False);"
```
