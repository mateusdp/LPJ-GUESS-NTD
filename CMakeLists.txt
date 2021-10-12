################################################################################
# CMake configuration file for building LPJ-GUESS
#
# To build LPJ-GUESS with this build system cmake needs to be installed.
# If it's not installed it can be downloaded for free from www.cmake.org.
#

cmake_minimum_required(VERSION 2.8.12.2...3.5.2)
project(guess)

# Compiler flags for building with Microsoft Visual C++
if (MSVC)
  # Disable warnings about using secure functions like sprintf_s instead of
  # regular sprintf etc.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _CRT_SECURE_NO_WARNINGS /D _SCL_SECURE_NO_WARNINGS")
endif()

# The following are configuration variables for the build. 
# Rather than editing them here, you should either edit them
# in the graphical cmake tools (cmake-gui on Windows, ccmake
# on Unix) or pass them in as parameters to cmake.

# A variable controlling whether or not to include unit tests
# Unit tests are disabled in old VS 6.0, since CATCH doesn't compile in such 
# an old compiler.
if (NOT MSVC_VERSION EQUAL "1200")
  set(UNIT_TESTS "OFF" CACHE BOOL "Whether to include unit tests")
endif()

if (UNIX)
  # Setup the SYSTEM variable, currently only used to choose which 
  # submit.sh to generate (for submitting to job queue)


  # Figure out what value it should have initially, based on the
  # environment variable ARCH if it's set.
  if (NOT $ENV{ARCH} STREQUAL "")
	set(DEFAULT_SYSTEM $ENV{ARCH})
  else()
	set(DEFAULT_SYSTEM "")
  endif()

  set(SYSTEM ${DEFAULT_SYSTEM} CACHE STRING "System to build for (empty (=simba), gimle, aurora or multicore)")
endif (UNIX)

# Where to search for cmake modules
# (used by cmake when the include() command is used in a cmake file)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# NetCDF - look for libraries and include files, and use them if found.
# If this CMakeList script does not find the path to the NetCDF library on
# your system, then consult LPJ-GUESS/reference/guessdoc.doc.
# On easybuild systems, it may help to use the same line as for aurora here.
# i.e.: set(CMAKE_PREFIX_PATH "$ENV{EBROOTNETCDF}").
if ("$ENV{ARCH}" STREQUAL "aurora")
  set(CMAKE_PREFIX_PATH "$ENV{EBROOTNETCDF}")
endif()
find_package(NetCDF QUIET)

if (NETCDF_FOUND)
  include_directories(${NETCDF_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${NETCDF_LIBRARIES})
  add_definitions(-DHAVE_NETCDF)
else()
  message("find_package: not found: NETCDF_C_LIBRARY and NETCDF_INCLUDE_DIR. See LPJ-GUESS documentation." )
  message("NETCDF disabled." )
endif()

# MPI - used if found (not needed on Windows)
if (NOT CMAKE_HOST_WIN32)
  find_package(MPI QUIET)
endif()

# These are deprecated according to documentation in the FindMPI module,
# but for some reason not hidden. Let's not display them for the typical
# LPJ-GUESS user who hasn't got MPI installed.
mark_as_advanced(MPI_LIBRARY MPI_EXTRA_LIBRARY)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  set(LIBS ${LIBS} ${MPI_LIBRARIES})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
  add_definitions(-DHAVE_MPI)
endif()

# Where the compiler should search for header files
include_directories(${guess_SOURCE_DIR}/framework ${guess_SOURCE_DIR}/libraries/gutil ${guess_SOURCE_DIR}/libraries/plib ${guess_SOURCE_DIR}/libraries/guessnc ${guess_SOURCE_DIR}/modules ${guess_SOURCE_DIR}/cru/guessio)

# The following directories contain source code and 
# additional CMakeLists.txt files
add_subdirectory(framework)
add_subdirectory(modules)
add_subdirectory(cru)
add_subdirectory(libraries)
if (WIN32)
  add_subdirectory(windows_version)
endif()

if (UNIT_TESTS)
  add_subdirectory(tests)
endif()

# Add the command line program's target
if (WIN32)
  # Let the exe be called guesscmd so it doesn't collide with the dll target
  set(guess_command_name "guesscmd")
else()
  # On Unix we don't have the dll so the command line binary can be called guess
  set(guess_command_name "guess")
endif()

# Specify the executable to build, and which sources to build it from
add_executable(${guess_command_name} ${guess_sources} command_line_version/main.cpp)

# Rule for building the unit test binary
if (UNIT_TESTS)
  add_executable(runtests ${guess_sources} ${test_sources})
  target_link_libraries(runtests ${LIBS})
endif()

# Specify libraries to link to the executable
target_link_libraries(${guess_command_name} ${LIBS})

if (WIN32)
  # Create guess.dll (used with the graphical Windows shell)
  add_library(guess SHARED ${guess_sources} windows_version/dllmain.cpp)

  # Specify libraries to link to the dll
  target_link_libraries(guess ${LIBS})
endif (WIN32)

# The custom build rule for generating the submit script from a template.
# The submit script is generated each time the command line binary is built.
if (UNIX)
  # Choose template file depending on system
  if (SYSTEM STREQUAL "")
    set(SYSTEM "simba")
  endif(SYSTEM STREQUAL "")
  set(submit_template ${guess_SOURCE_DIR}/parallel_version/${SYSTEM}.tmpl)

  # A sed command which exchanges each occurrence of BINARY and DATE in
  # the template with appropriate path to the compiled binary and the
  # current date.
  add_custom_command(TARGET ${guess_command_name}
                     POST_BUILD
             COMMAND sed 's|BINARY|'${CMAKE_CURRENT_BINARY_DIR}/${guess_command_name}'|g\;s|DATE|'\"`date`\"'|g' ${submit_template} > ${CMAKE_BINARY_DIR}/submit.sh
			 COMMAND chmod +x ${CMAKE_BINARY_DIR}/submit.sh
             COMMENT "Generating submit script")
endif (UNIX)

# Rule for running unit tests automatically
if (UNIT_TESTS)
  add_custom_command(TARGET runtests 
    POST_BUILD 
    COMMAND runtests 
    COMMENT "Running tests")
endif()

if (UNIX)
   # pgCC 6 doesn't seem to recognize -rdynamic, so remove it
   # (we shouldn't need it anyway)
   # It seems the CMake developers have fixed this in newer versions
   # (sometime after 2.8)
   SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
   SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
endif(UNIX)

if (UNIX)
  # Set default build type to Release on Unix
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
  endif(NOT CMAKE_BUILD_TYPE)
endif(UNIX)
