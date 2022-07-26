#!/bin/bash
# set -x
#
# Driver for QuICC benchmarking
#


# default params
quicc_root="../QuICC"
launch_script="echo"
max_tasks_per_node="36" # mc node
min_tasks_per_node="1" # number of sockets
max_threads_per_task="1"    # no threads
steps="100"
time="00:15:00"

# version
version="0.0.1"
version_date="17-05-2022"
version_string="profiling_driver $version ($version_date)"

# show usage and handle arguments
function showusage {
    echo "\
Driver for benchmarking.

USAGE:
    profiling_driver [OPTIONS].

OPTIONS:
    -p PATH, --path PATH
            QuICC root directory.
    -l SCRIPT, --launch SCRIPT
            Script to run simulation.
    -t HH:MM:SS, --time HH:MM:SS
            Time to run simulation.
    -M TASKS, --max_tasks_per_node TASKS
            Max tasks per node.
    -m TASKS, --min_tasks_per_node TASKS
            Min tasks per node.
    -mt THREADS, --max_threads_per_task THREADS
            Max thread per task.
    -s STEPS, --steps STEPS
            Number of time steps.
    -d, --dry-run
    -h, --help
            Display this usage information.
    -V, --version
            Show version and exit.
"
}

# parse args
while [[ "$1" != "" ]]; do
    case $1 in
        "")
            shift
            ;;
        -p|--path)
            quicc_root=$2
            shift
            shift
            ;;
        -l|--launch)
            launch_script=$2
            shift
            shift
            ;;
        -t|--time)
            time=$2
            shift
            shift
            ;;
        -M|--max_tasks_per_node)
            max_tasks_per_node=$2
            shift
            shift
            ;;
        -m|--min_tasks_per_node)
            min_tasks_per_node=$2
            shift
            shift
            ;;
        -mt|--max_threads_per_task)
            max_threads_per_task=$2
            shift
            shift
            ;;
        -s|--steps)
            steps=$2
            shift
            shift
            ;;
        -d|--dry-run)
            echo "debug mode"
            debug=echo
            shift
            ;;
        -V|--version)
            echo "$version_string"
            exit 0
            ;;
        -h|-\?|--help)
            showusage
            exit 0
            ;;
        *)
            showusage
            exit 1
            ;;
    esac

done

# QuICC versions
# needs to match build-${quicc_version}
declare -a quicc_version=(
    # "gcc"
    "gcc-serial-openmp"
    )

declare -a problem_size=(
    "47;95;95"
    "63;127;127"
    )

declare -a problem_tasks=(
    "4"
    "8"
    "16"
    "32"
    "36"
    # "72"
    # "144"
    )

# QuICC model
quicc_model="BoussinesqSphereDynamo"

declare -a input_files=(
    "parameters.cfg" \
    "state_initial.hdf5"
    )

# loop over versions
for qv in ${quicc_version[@]}; do
    echo QuICC version: ${qv}

    # loop over problem size
    for size in ${problem_size[@]}; do
        # unpack size
        IFS=";" read -r -a arr <<< "${size}"
        nr="${arr[0]}"
        nl="${arr[1]}"
        nm="${arr[2]}"

        echo "nr: ${nr}, nl: ${nl}, nm: ${nm}"

        # loop over tasks
        for tasks in ${problem_tasks[@]}; do
            echo tasks: $tasks

            # get nodes and tasks per node
            nodes=$(($tasks / $max_tasks_per_node / $max_threads_per_task))
            nodes=$(($nodes > 0 ? $nodes : 1))
            echo nodes: $nodes

            ntasks_per_node=$(($tasks / $nodes / $max_threads_per_task))
            floor_min=$(($min_tasks_per_node > 1 ? $min_tasks_per_node - 1 : 0))
            ntasks_per_node=$(($ntasks_per_node > $floor_min ? $ntasks_per_node : $min_tasks_per_node))
            echo ntasks_per_node: $ntasks_per_node

            threads_per_task=$(($tasks / $nodes / $ntasks_per_node))
            echo threads_per_task: $threads_per_task


            # create and move to folder
            experiment_folder="${qv}/${quicc_model}/size_${nr}_${nl}_${nm}/tasks_${tasks}"
            if [ -d "${experiment_folder}" ]; then
                echo ${experiment_folder} already exists, skipping
                continue
            fi
            ${debug} mkdir -p -- "${experiment_folder}" && ${debug} pushd "$_"

            # copy input files
            ref_path="${quicc_root}/build-${qv}/Models/${quicc_model}/TestSuite/Benchmarks/_refdata/Explicit"
            for file in "${input_files[@]}"; do
                ${debug} cp "${ref_path}/${file}" .
            done

            # adjust input process
            ${debug} sed -i "s/<cpus>[0-9]\+<\/cpus>/<cpus>${tasks}<\/cpus>/" parameters.cfg
            # adjust input size
            ${debug} sed -i "s/<dim1D>[0-9]\+<\/dim1D>/<dim1D>${nr}<\/dim1D>/" parameters.cfg
            ${debug} sed -i "s/<dim2D>[0-9]\+<\/dim2D>/<dim2D>${nl}<\/dim2D>/" parameters.cfg
            ${debug} sed -i "s/<dim3D>[0-9]\+<\/dim3D>/<dim3D>${nm}<\/dim3D>/" parameters.cfg
            # adjust number of steps
            ${debug} sed -i "s/<sim>-?[0-9]\+<\/sim>/<sim>-${steps}<\/sim>/" parameters.cfg
            # disable output
            ${debug} sed -i "s/<ascii>[0-9]\+<\/ascii>/<ascii>0<\/ascii>/" parameters.cfg
            ${debug} sed -i "s/<hdf5>[0-9]\+<\/hdf5>/<hdf5>0<\/hdf5>/" parameters.cfg
            ${debug} sed -i "s/<enable>[0-9]\+<\/enable>/<enable>0<\/enable>/g" parameters.cfg

            # symbolic link of executable
            ${debug} ln -s "${quicc_root}/build-${qv}/bin/${quicc_model}ExplicitModel"


            # launch
            ${debug} sbatch --nodes=${nodes} \
                --ntasks-per-node=${ntasks_per_node} \
                --cpus-per-task=${threads_per_task} \
                --job-name=${qv}_${quicc_model} \
                --time=${time} \
                 ${launch_script}

            # go back to base folder
            ${debug} popd

        done # tasks loop

    done # size loop

done # version loop


# to do, sweep problem size
