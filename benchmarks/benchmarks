#!/bin/bash
#
# benchmarks - an automatic benchmarking suite for LPJ-GUESS.
#
# For information on how to use, see the readme.txt, or run this
# script without arguments.
#

cmake_command="cmake -D CMAKE_BUILD_TYPE=Release"
make_command="make -j2"
bbox_command="bbox"
pp_walltime=5:00:00

# Prints a message to stderr
function error {
    echo "$1" >&2
}

# A function which attempts to find out whether PBS is installed
function has_pbs {
    # Assume PBS is installed correctly if the qsub command is found
    type -P qsub &>/dev/null
}

# A function which attempts to find out whether SLURM is installed
function has_slurm {
    # Assume SLURM is installed correctly if the sbatch command is found
    type -P sbatch &>/dev/null
}

# Informs the user of how to use the script, and then exits.
function print_usage_and_exit {
    error "Usage: benchmarks [ -s ] [ -c ] [ -i TESTS | -x TESTS ] <workdir>"
    error
    error "Options"
    error "======="
    error "-s   Single machine mode, don't use PBS"
    error "-c   Compile only, don't run the tests"
    error "-p   Post processing only"
    error "-r   Generate printer friendly report"
    error "     (does nothing else and is run as a regular process,"
    error "      in other words, not submitted as a batch job)"
    error "-i   Include only the specified tests"
    error "-x   Exclude the specified tests"
    error
    error "Examples"
    error "========"
    error "1. Run all benchmarks:"
    error "benchmarks /scratch/sally/benchmarks"
    error "2. Run some benchmarks:"
    error "benchmarks -i \"europe euroflux\" /scratch/sally/benchmarks"
    error "3. Run all except some benchmarks:"
    error "benchmarks -x \"global europe\" /scratch/sally/benchmarks"
    exit 1
}

# A help function to check if something is included in a space separated list
# (used to check if a test name is included in a list of tests to include or
# exclude)
# Parameters:
# $1  A list of names
# $2  A name to search for
function includes {
    for name in $1; do
	if [[ $name == $2 ]]; then
	    return 0
	fi
    done
    return 1
}

single_machine=0
compile_only=0
pp_only=0
printer_friendly=0

# Handle the -s, -c, -i and -x options
while getopts ":scpri:x:" opt; do
    case $opt in
	i  ) inclusive=$OPTARG
	     inclusive_flag=1 ;;
	x  ) exclusive=$OPTARG
	     exclusive_flag=1 ;;
	s  ) single_machine=1 ;;
	c  ) compile_only=1 ;;
	p  ) pp_only=1 ;;
	r  ) printer_friendly=1 ;;
	\? ) print_usage_and_exit ;;
    esac
done

if (( ${#inclusive} > 0 && ${#exclusive} > 0 )); then
    error "Using both -i and -x is not allowed!"
    exit 1
fi

if (( $single_machine == 0 )); then
    if ! has_pbs && ! has_slurm ; then
	error "It doesn't seem like PBS or SLURM is installed. You can run the"
	error "benchmarks on the local machine with the -s option."
	exit 1
    fi
else
	# Use the multicore submit script if running on a local machine
	cmake_command="$cmake_command -D SYSTEM=multicore"
fi

if (( $compile_only == 1 && $pp_only == 1 )); then
    error "-c and -p can't be used at the same time!"
    exit 1
fi

if (( $printer_friendly == 1 && ( $compile_only == 1 || $pp_only == 1) )); then
    error "-r can't be combined with -c or -p!"
    exit 1
fi

# Throw away the options getopts handled
shift $(($OPTIND -1))

# Make sure the user supplied a working directory
if [ $# != 1 ]; then
    print_usage_and_exit
fi

work_dir=$1

# The build function is used to build one of the versions of the model.
# Parameters:
# $1  The directory where the results of the build should be placed
# $2  (optional) Path to additional source code
function build {
    local build_path=$1
    local additional_source=$2

    echo "Building in $build_path"

    rm -rf $build_path
    mkdir -p $build_path/output

    source_dirs="framework modules cru libraries command_line_version parallel_version cmake"
    # Copy with tar so we can exclude svn directories
    (cd ..; tar cf - --exclude=".svn" $source_dirs) | (cd $build_path; tar xf -)

    cp ../CMakeLists.txt $build_path

    if [ -n "$additional_source" ]; then
	# Copy with tar so we can exclude svn directories
	(cd $additional_source; tar cf - --exclude=".svn" *) | (cd $build_path; tar xf -)
    fi

    pushd $build_path/output &> /dev/null

    cmake_output=$(mktemp)
    if ! $cmake_command .. &> $cmake_output ; then
	error "Failed to configure build " $build_path "!"
	cat $cmake_output
	rm $cmake_output
	popd &> /dev/null
	return 1
    fi
    rm $cmake_output

    make_output=$(mktemp)
    if ! $make_command &> $make_output ; then
	error "Failed to build " $build_path "!"
	cat $make_output
	rm $make_output
	popd &> /dev/null
	return 1
    fi
    rm $make_output

    popd &> /dev/null
    return 0
}

# This function prepares the workspace for one of the benchmarks
# by copying in necessary files and creating the submit script.
# Parameters:
# $1  The name of the benchmark
# $2  Path to the build directory
function setup_workspace {
    local test_name=$1
    local build_path=$2

    local workspace=$work_dir/$test_name

    rm -rf $workspace
    mkdir $workspace

    # Get the submit script
	cp $build_path/output/submit.sh $workspace

    # Copy the config files to the workspace.
    # The ins file is excluded to avoid confusion. It's the ins
    # file in the benchmark directory that's used. If we copy
    # the ins file to the run directory it looks like the ins
    # file in the benchmark directory can be modified once the
    # job has been submitted to queue.
    ls $test_name/config/* | grep -v guess.ins | xargs -I '{}' cp '{}' $workspace

    # Limit the gridlist if there is a boundingbox.txt
    if [ -f $test_name/config/boundingbox.txt ]; then
	bbox=$(cat $test_name/config/boundingbox.txt)
	$bbox_command -f $test_name/config/gridlist.txt \
	    -o $workspace/gridlist.txt $bbox
    fi

	# Create a symbolic link to the common directory
	ln -s -T $(readlink -f common) $workspace/common
}

# This function submits the post processing job.
# Parameters:
# $1  The name of the benchmark
# $2  Job-id of the actual benchmark job (for dependency)
#     Don't supply if running only the post processing
function submit_pp {
    local test_name=$1
    local workspace=$work_dir/$test_name
    local bm_root=$PWD

    # Construct argument for qsub/sbatch for dependency if there is one
    if [ -n "$2" ]; then
	if has_pbs ; then
	    dependency="-W depend=afterok:$jobid"
	else
	    dependency="--dependency=afterok:$jobid"
	fi
    else
	dependency=""
    fi

    pushd $workspace &> /dev/null
    # Create the post processing job

    # The part of the job script which is the same for PBS and SLURM
    common_script=$(mktemp)
    cat <<EOF > $common_script
if [ -f $bm_root/$test_name/postprocess.sh ]; then
    PATH=$PATH $bm_root/runpostprocessing $bm_root/$test_name/postprocess.sh $test_name
fi
ls | grep "^run[0-9]*$" | xargs rm -rf
cd $bm_root
echo '$test_name' >> progress.txt
EOF

    if has_pbs ; then
	# Submit as a PBS job
	qsub -l walltime=$pp_walltime $dependency <<EOF &> /dev/null
#!/bin/sh
#PBS -N bm_${test_name}_pp
cd $workspace
$(cat $common_script)
EOF
    else
	# Submit as a SLURM job
	sbatch -n 1 --time=$pp_walltime $dependency <<EOF &> /dev/null
#!/bin/sh
#SBATCH -J bm_${test_name}_pp
$(cat $common_script)
EOF
    fi

    rm $common_script

    popd &> /dev/null    
}

# This function submits a benchmark as a PBS job
# Parameters:
# $1  The name of the benchmark
function submit_job {
    local test_name=$1
    local workspace=$work_dir/$test_name
    local bm_root=$PWD
    local submit_vars=$PWD/$test_name/submit_vars.sh
    local insfile=$bm_root/$test_name/config/guess.ins

    pushd $workspace &> /dev/null
    if [ -f $submit_vars ]; then
	jobid=$(./submit.sh -i $insfile -n "bm_$test_name" -s $submit_vars 2> /dev/null)
    else
	jobid=$(./submit.sh -i $insfile -n "bm_$test_name" 2> /dev/null)
    fi

    local submit_exit_status=$?

    popd &> /dev/null

    if (($submit_exit_status != 0)); then
	error "Failed to submit $test_name! Try manually for more information."
    else
	submit_pp $test_name $jobid
    fi
}

# This function runs the post processing as a regular process
# (not using PBS)
# Parameters:
# $1 The name of the benchmark
function run_single_machine_pp {
    local test_name=$1
    local workspace=$work_dir/$test_name
    local bm_root=$PWD

    pushd $workspace &> /dev/null
    if [ -f $bm_root/$test_name/postprocess.sh ]; then
	nice $bm_root/runpostprocessing $bm_root/$test_name/postprocess.sh $test_name
    fi
    ls | grep "^run[0-9]*$" | xargs rm -rf
    popd &> /dev/null

    echo $test_name >> $bm_root/progress.txt
}

# This function runs a benchmark as a regular process
# (not using PBS)
# Parameters:
# $1 The name of the benchmark
function run_single_machine_model {
    local test_name=$1
    local workspace=$work_dir/$test_name
    local bm_root=$PWD
    local submit_vars=$PWD/$test_name/submit_vars.sh
    local insfile=$bm_root/$test_name/config/guess.ins

    pushd $workspace &> /dev/null
    if [ -f $submit_vars ]; then
		./submit.sh -i $insfile -s $submit_vars
    else
		./submit.sh -i $insfile
    fi

    popd &> /dev/null

    run_single_machine_pp $test_name
}

# This function takes care of everything that needs to be done 
# for a benchmark: building, setting up the workspace,
# submitting the model job and submitting the post processing job.
# Parameters:
# $1  The test name
function run_benchmark {
    local test_name=$1
    local build_path=$work_dir/standard_build

    echo $test_name"..."

    if [ -f $test_name/extra_source.txt ]; then
	build_path=$work_dir/${test_name}_build
	echo "Building $test_name model..."
	if ! build $build_path extra_source/$(cat $test_name/extra_source.txt); then
	    echo "Not running ${test_name}..."
	    return
	fi
    fi

    if (( $compile_only == 1 )); then
	return
    fi

    setup_workspace $test_name $build_path

    if (( $single_machine == 1 )); then
	run_single_machine_model $test_name
    else
	submit_job $test_name
    fi
}

# This function runs just the post processing for one benchmark.
# Parameters:
# $1 The test name
function run_pp {
    local test_name=$1
    local workspace=$work_dir/$test_name

    if [ -d $workspace ]; then
	if (( $single_machine == 1 )); then
	    run_single_machine_pp $test_name
	else
	    submit_pp $test_name
	fi
    else
	echo "$workspace doesn't exist, skipping!"
    fi
}

# This function prints out interesting information about this benchmarks
# run to a text file in the root of the working directory.
# Parameters:
# $1 Absolute path to the working directory
function meta_data {
    file=$1/info.txt
    file2=$1/info.local_modifications.txt

    echo "Started at: "`date` > $file
    echo "Hostname:   "`hostname` >> $file
    echo "User:       "`whoami` >> $file

    echo >> $file

    # Go up to parent directory and print out version control information
    pushd .. &> /dev/null

    if svn info &> /dev/null ; then
	# Strip host name from svn URL
	svn info | grep "^URL:" | sed 's|URL: svn://.*/svn/|VCS Path: |' >> $file
	svn info | grep "^Revision:" >> $file
	svn info | grep "^Last Changed Date:" >> $file
	
	changes=$(svn st -q --non-interactive | grep -v "Performing status on external item" | grep -v "^$")
	changesdiff=$(svn diff --non-interactive)

	if [ -n "$changes" ]; then
	    echo >> $file
	    echo "Local modifications:" >> $file
	    echo "$changes" >> $file

	    # Print local WC path and svn diff
	    echo "Local WC path:" > $file2
	    pwd -P >> $file2
	    echo >> $file2
	    echo "$changesdiff" >> $file2
	fi
    else
	echo "Version information missing" >> $file
    fi

    popd &> /dev/null
}

# Ok, let's start!

# Create the progress.txt file, overwriting if it existed
echo "The following benchmarks are done:" > progress.txt

# If necessary, create the working directory where all the benchmarks will be
# run from.
mkdir -p $work_dir

# Convert the work_dir variable to an absolute path if it isn't already
work_dir=`cd $1; pwd`

if (( $printer_friendly == 1 )); then
    ./printerfriendly $work_dir
    exit 0
fi

# Generate meta data about this run 
# (such as time stamp and version control information).
if (( $pp_only == 0 )); then
    meta_data $work_dir
fi

# Build the model which will be used by benchmarks that don't have their own
# special code. No need to build it if we're only doing post processing
# however.
if (( $pp_only == 0 )); then
    echo "Building standard model..."
    if ! build $work_dir/standard_build; then
	echo "Not running any benchmarks!"
	exit 1
    fi
fi

# Find all sub directories
subdirs=($(find . -mindepth 1 -maxdepth 1 -type d))

# Figure out which tests to run
num_subdirs=${#subdirs[@]}
for ((i=0; i < num_subdirs; i++));
do
  subdir=${subdirs[i]}
  # It needs to have a config subdir to be considered
  if [ -d $subdir/config ]; then
      # Remove "./" to convert from path to test name
      test_name=$(basename $subdir)

      # Include the test according to -x or -i flag
      if [ -n "$inclusive_flag" ]; then
	  if includes "$inclusive" $test_name; then
	      test_names[i]=$test_name
	  fi
      elif [ -n "$exclusive_flag" ]; then
	  if ! includes "$exclusive" $test_name; then
	      test_names[i]=$test_name
	  fi
      else
	  # No flag, include unconditionally
	  test_names[i]=$test_name
      fi
  fi
done

# Call run_benchmark for each test
for test_name in "${test_names[@]}"
do
  if (( $pp_only == 1 )); then
      run_pp $test_name
  else
      run_benchmark $test_name
  fi
done

if (( $compile_only == 1 )); then
    echo "Done!"
    exit 0
fi

if (( $single_machine == 0)); then
    echo "All jobs submitted!"
    echo "See the file progress.txt for information of completed benchmarks."
fi
