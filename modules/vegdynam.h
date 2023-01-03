///////////////////////////////////////////////////////////////////////////////////////
/// \file vegdynam.h
/// \brief Vegetation dynamics and disturbance
///
/// \author Ben Smith
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_VEGDYNAM_H
#define LPJ_GUESS_VEGDYNAM_H

#include "guess.h"

void vegetation_dynamics(Stand& stand,Patch& patch);
void planting(Patch& patch);

#endif // LPJ_GUESS_VEGDYNAM_H
