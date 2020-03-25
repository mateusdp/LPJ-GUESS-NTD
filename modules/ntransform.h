///////////////////////////////////////////////////////////////////////////////////////
/// \file ntransform.h
/// \brief Framework header file, LPJ-GUESS Combined Modular Framework
///
/// This header file contains:
///  (1) definitions of all main classes used by the framework and modules. Modules may
///      require classes to contain certain member variables and functions (see module
///      source files for details).
///  (2) other type, constant and function definitions to be accessible throughout the
///      model code.
///  (3) a forward declaration of the framework function if this is not the main
///      function.
///
/// \author Xu-Ri and modified for LPJ-guess by Peter Eliasson, David Wårlind and Stefan Olin.
/// $Date: 2013-10-14 14:12:00 +0100 (Mon, 10 Sep 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Nitrogen transformation processes in Soil
//
//
// Header file name:      ntransform.h
// Source code file name: ntransform.cpp
// Written by:            Stefan Olin, adopted from Xu-Ri 2004-06-10.
// Version dated:         2019.
//
// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_NTRANSFORM_H
#define LPJ_GUESS_NTRANSFORM_H

#include "guess.h"

void ntransform(Patch& patch, Climate& climate);

#endif // LPJ_GUESS_NTRANSFORM_H
