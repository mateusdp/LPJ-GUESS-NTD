////////////////////////////////////////////////////////////////////////
/// \file soilmethane.h
/// \brief Declarations used in methane code
/// The class Soil and its member functions and variables are declared in guess.h, 
/// while its member functions are implemented in soil.cpp and in soilmethane.cpp.
/// 
/// \author Paul Miller
/// $Date: 2019-04-04 13:44:06 +0200 (Thu, 04 Apr 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_METHANE_H
#define LPJ_GUESS_METHANE_H

#include "config.h"
#include "guess.h"
#include "soil.h"

/// Methane dynamics for this patch today
/** Calculates methane fluxes from each soil layer from various emission pathways
 *  Called each day for wetland stands.
 */
void methane_dynamics(Patch& patch);

#endif // !LPJ_GUESS_METHANE_H
