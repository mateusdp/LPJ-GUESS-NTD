///////////////////////////////////////////////////////////////////////////////////////
/// \file config.h
/// \brief Configuration header file, accessible throughout the model code.
///
/// This header file contains preprocessor definitions accessible throughout the
/// model code. Each implementation file (.cpp file) in LPJ-GUESS should include
/// this header before anything else. Definitions should be placed here if they
/// are of interest when configuring how to build LPJ-GUESS, for instance
/// definitions for enabling/disabling a module, or changing a behaviour which
/// for some reason isn't configurable from the instruction file.
///
/// This file may also contain non-model related code for working around platform
/// specific issues, such as non-standard conforming compilers.
///
/// \author Joe Siltberg
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CONFIG_H
#define LPJ_GUESS_CONFIG_H

// Compiler specific checks, for instance for disabling specific warnings

// All versions of Microsoft's compiler
#ifdef _MSC_VER

// 'this' : used in base member initializer list
#pragma warning (disable: 4355)
// long name
#pragma warning (disable: 4786)

#endif

// min and max functions for MS Visual C++ 6.0
#if defined(_MSC_VER) && _MSC_VER == 1200
namespace std {
template <class T> inline T max(const T& a, const T& b) {
    return (a > b) ? a : b;
}

template <class T> inline T min(const T& a, const T& b) {
    return (a < b) ? a : b;
}
}
#else
// modern compilers define min and max in <algorithm>
#include <algorithm>
#endif

using std::min;
using std::max;

// platform independent function for changing and making working directory
// we'll call our new functions change_directory and make_directory
// make_directory is non recursive.
#ifdef _MSC_VER
// The Microsoft way
#include <direct.h>
#define change_directory _chdir
static void make_directory(const char *dir) {
	_mkdir(dir);
}
#else
// The POSIX way
#include <unistd.h>
#include <sys/stat.h>
#define change_directory chdir
static void make_directory(const char *dir) {
	mkdir(dir, 0755);
}
#endif

#endif // LPJ_GUESS_CONFIG_H
