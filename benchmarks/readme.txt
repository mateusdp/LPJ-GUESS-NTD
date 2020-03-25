

                   Automatic benchmarking suite for LPJ-GUESS


This directory contains benchmarks and tools for running them. To run the
benchmarks, the 'benchmarks' script in this directory is used. Each benchmark
has its own directory where configuration can be done.

Running the benchmarks
======================
To run the benchmarks you need a Unix system. You can run the tests as regular
processes, in which case most Unix systems should work. By default however the
tests are submitted to a PBS batch queue. This has been tested on Simba 
(ENES own cluster), but it should be possible to adapt the scripts for any
system with PBS installed.

Dependencies
------------
The post processing scripts which are run after each finished simulation 
currently depend on the following software:

 - guess_utilities
 - gnuplot
 - ImageMagick (only needed if printer friendly reports are to be generated)

They are only needed for the post processing. If missing, the simulations will
still run but the post processing will fail. 

guess_utilities can be found at the LPJ-GUESS download page, or in the version 
control repository. The others are standard Unix software packages, and are
often installed by default.

Basic usage
-----------
The benchmarks script should be run from this directory. When starting it you
need to supply a working directory from where the benchmarks will run. On 
Simba this directory should be somewhere under scratch.

The easiest way to run the benchmarks is like this:

./benchmarks /scratch/sally/benchmarks

This will run all available benchmarks, placing all files under
/scratch/sally/benchmarks. The produced output files from for instance the
global benchmark will be under /scratch/sally/benchmarks/global.

For each benchmark two PBS jobs are submitted, one for the actual model, and
one for post processing of the output. You can use PBS tools like qstat to
see when your jobs are finished, or you can have a look in the file 
progress.txt in this directory.

The script allows you to specify exactly which benchmarks to run, and has
an option for running without PBS. For more information, run the script
without arguments.

Adding or configuring benchmarks
--------------------------------
Each benchmark is defined by the contents of a directory. A minimal benchmark
looks like this:

mybenchmark
|-- config
    |-- gridlist.txt
    `-- guess.ins

The benchmark will be named by the directory, in this case 'mybenchmark'.
Each benchmark needs to have a config directory with the files gridlist.txt
and guess.ins which are used to configure the model. In addition, you can
add optional files, an example with all possible additional configuration
may look like this:

mybenchmark
|-- config
|   |-- gridlist.txt
|   `-- guess.ins
|-- extra_source.txt
|-- postprocess.sh
`-- submit_vars.sh

postprocess.sh is a bash script file containing post processing commands 
specific for this benchmark. 

If there is a file named extra_source.txt, it should include the name of a
directory under benchmarks/extra_source. The files in that directory will be 
included in the compilation. Additional source should be avoided when possible 
since it can easily become a maintenance problem. The source files in this 
directory must follow exactly the same structure as in the regular source code. 
For instance if you want to replace modules/driver.cpp your replacement must
be named driver.cpp and be placed in a directory named modules.

If there is a file named submit_vars.sh it will be used to configure the PBS 
jobs. The file can define the variables NPROCESS, WALLTIME, INSFILE, 
GRIDLIST and OUTFILES, for instance like this:

NPROCESS=20
WALLTIME=12:00:00

Where:
NPROCESS = number of processes in parallel job
WALLTIME = maximum wall (real) time for job hh:mm:ss
INSFILE  = path to ins file from run directory
GRIDLIST = path to gridlist file from run directory
OUTFILES = list of LPJ-GUESS output files in single quotes,
           and separated by spaces (filenames only, including
           extension, no directory.)

If a variable isn't defined it will be given a default value.

Post processing
---------------
Post processing commands that should be performed for a benchmark should be
placed in the file postprocess.sh in the benchmark's directory. In order to
avoid duplication, this script can call other scripts from the "postprocess"
directory (in the "benchmarks" directory). The scripts under postprocess
contain commands that can be useful in several of the benchmarks.

Common files
------------
Files used by several benchmarks, for instance legends used by the post 
processing etc., can be placed in the "common" directory. This directory
is accessible from each benchmark and post processing script through
a symbolic link in the benchmarks working directory.


Joe Siltberg
joe.siltberg@nateko.lu.se
2011-06-10
