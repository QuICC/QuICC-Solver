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
reframe -C $QUICC_ROOT/ci/reframe/settings.py -c $QUICC_ROOT/ci/reframe/quicc_library.py -r -S quicc_root=$QUICC_ROOT --performance-report --force-local
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
reframe -c $QUICC_ROOT/ci/reframe/quicc_library.py -r -S quicc_root=$QUICC_ROOT --performance-report --force-local
```
