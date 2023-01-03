///////////////////////////////////////////////////////////////////////////////////////
/// \file main.cpp
/// \brief Main module for command line version of LPJ-GUESS
///
/// \author Ben Smith
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"
#include "framework.h"
#include "commandlinearguments.h"
#include "parallel.h"
#include <stdlib.h>

#ifdef __unix__
// includes needed for umask()
#include <sys/types.h>
#include <sys/stat.h>
#endif

///////////////////////////////////////////////////////////////////////////////////////
// LOG FILE
// The name of the log file to which output from all dprintf and fail calls is sent is
// set here

xtring file_log="guess.log";



///////////////////////////////////////////////////////////////////////////////////////
// MAIN
// This is the function called when the executable is run

int main(int argc,char* argv[]) {
	// Initialize parallel communication if available.
	// This needs to be done before command line parsing since some MPI
	// implementations put their own arguments in our command line
	// (which should then be removed after the MPI initialization).
	GuessParallel::init(argc, argv);

#ifdef __unix__

	// On unix systems we set the umask (which controls the file permissions
	// of files we create). Some MPI implementations set a restrictive umask
	// which would mean that other users can't read the files generated by
	// LPJ-GUESS.

	umask(S_IWGRP | S_IWOTH); // only disable write access for group and others

#else

	// Maximizing capacity of number of open files
	_setmaxstdio(2048);

#endif

	// Parse command line arguments
	CommandLineArguments args(argc, argv);

	// Change working directory according to rank if parallel run
	if (args.get_parallel()) {
		xtring path;
		path.printf("./run%d", GuessParallel::get_rank()+1);

		if (change_directory(path) != 0) {
			fprintf(stderr, "Failed to change to run directory\n");
			return EXIT_FAILURE;
		}
	}

	// Set our shell for the model to communicate with the world
	set_shell(new CommandLineShell(file_log));

	if (args.get_help()) {
		printhelp();
		return EXIT_SUCCESS;
	}

	// Call the framework
	framework(args);

	// Say goodbye
	dprintf("\nFinished\n");

	return EXIT_SUCCESS;
}
