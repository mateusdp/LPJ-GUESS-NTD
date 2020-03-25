///////////////////////////////////////////////////////////////////////////////////////
/// \file somdynam.h
/// \brief Soil organic matter dynamics
///
/// \author Ben Smith
/// $Date: 2017-04-24 19:33:38 +0200 (Mo, 24. Apr 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SOMDYNAM_H
#define LPJ_GUESS_SOMDYNAM_H

#include "guess.h"

void som_dynamics(Patch& patch);

#endif // LPJ_GUESS_SOMDYNAM_H
