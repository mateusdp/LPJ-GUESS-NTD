                        LPJ-GUESS Version 4.1
                        =====================

                        PLEASE READ CAREFULLY
 
This directory contains source code files and other files necessary to
run a *demonstation version* of LPJ-GUESS:
     - in population mode (i.e. as LPJ-DGVM) or cohort mode (i.e. as GUESS)
     - on the Unix or Windows operating system
     - in Windows, as either a command-line executable or within the
       LPJ-GUESS Windows Shell

Note that the concept of "plug and play" does not apply to LPJ-GUESS.
It is the responsibility of each individual user to:
     - provide a C++ compiler. For installation under Windows a software
       development "platform" is also strongly recommended. Documentation
       is available for Microsoft Visual C++.
     - optain the appropriate climatic, atmospheric CO2 and soil data
       to drive the model for their particular study
     - provide source code to read their particular input data and
       produce their particular output
     - produce an instruction script (ins file) with run-time parameters
       (data file names, model configuration settings, PFT parameter values
       etc) for the model runs
Some technical advice regarding these matters can be found in the draft
documentation (reference/guessdoc.pdf). Additional help is available as
commenting in model source code and header files.

The model is compiled with the CMake build system. If it is not available
on your system you can download it for free from www.cmake.org. You can
find further information on how to compile the model with CMake in the
draft documentation (reference/guessdoc.pdf).


IMPORTANT: FTP-TRANSFERS FROM UNIX TO WINDOWS:
The following files should be transferred in "ASCII" mode:
     - Files with the extensions .cpp, .h, .grd, .txt, .dat and .ins
     - Makefile
The following files should be transferred in "binary" mode:
     - Files with the extensions .lib, .exe, .pdf, .bin and .gz
     

Structure of this directory:

./reference
     - Current draft of LPJ-GUESS technical manual

./framework
     - Model framework source code and header files (for explanation see
       technical manual)

./modules
     - Source code and header files for model modules (except "Main"
       module). An input module (demoinput.*) is provided for demonstration
       purposes, compatible with the input data supplied in
       directory ../data. There's also an input module for NetCDF data
       (must comply with the CF metadata standard). Most users will be able
       to write an input module customised to their own study by modifying 
       the demonstration version supplied. A standard output module is also
       provided (commonoutput.*), with many standard outputs.
       Further explanations in the technical manual.

./libraries
     - Source code and header files for the custom libraries gutil, plib, 
       and guessnc required by LPJ-GUESS.

./command_line_version
     - Files required specifically to install LPJ-GUESS as a command-line
       executable on Unix or Windows.

./windows_version
     - Files required specifically to install LPJ-GUESS as a dynamic link
       library (DLL) to run under the LPJ-GUESS Windows Shell. This is the
       recommended configuration for running the model under Windows.

./data
     - Input data files for the demonstration version of LPJ-GUESS.
       Subdirectories contain:
       ./env
            - historical climate data, soil data. See "readme" file.
       ./ins
            - sample instruction script (ins) files for global and european
              simulations, compatible with demonstration input. 
              See commenting in files.
       ./gridlist
            - sample grid cell coordinate list files compatible with
              demonstration input module and instruction script
              files under directory ../ins
      
./benchmarks
     - Files required for running the benchmarks (under linux only).
       Note that the postprocessing of the benchmarks requires 
	   the linux helper programs of "LPJ-GUESS utilities": its source 
	   code is available at the same www page as where this version of 
       LPJ-GUESS was downloaded from. The source code needs to be built 
	   on the system where you will run it. Do not move executables
	   after build.
              
./cru
     - Input module version for reading in CRU-NCEP historical climate 
       data for 1901-2015 in custom binary format used by LPJ-GUESS.
       The file name cru_TS30.cpp and the file name cru_TS30.cpp/.h 
       (and cognate namespace CRU_TS30) has been kept,  although the 
       forcing data is no longer CRU TS.
       
       NOTE: The associated CRU-NCEP binary files are not in the public domain. 
       Normally users should create their own "./cru" binary files, e.g. using  
       the Fast Archive project available the same level as where this version of 
       LPJ-GUESS was downloaded from.
       
Joe Siltberg 2014-06-27 (Version 3.0)
Johan Nord 2021-09-29 (Version 4.1)
