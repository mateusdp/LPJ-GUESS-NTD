///////////////////////////////////////////////////////////////////////////////////////
/// \file somdynam.h
/// \brief Soil organic matter dynamics
///
/// \author Ben Smith
/// $Date: 2021-04-22 18:36:50 +0200 (Thu, 22 Apr 2021) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SOMDYNAM_H
#define LPJ_GUESS_SOMDYNAM_H

#include "guess.h"

// computes the decomposition temperature modifier (in range 0-1)
double temperature_modifier(double temp);

// computes the decomposition moisture modifier (in range 0-1)
double moisture_modifier(double wfps);

void som_dynamics(Patch& patch, Climate& climate);

// computes the fraction of leaf and root that goes to metabolic litter (used by BLAZE)
double metabolic_litter_fraction(double lton);

// determine lignin to N ratio for leaf and root litter (used by BLAZE)
double lignin_to_n_ratio(double cmass_litter, double nmass_litter, double LIGCFRAC, double cton_avr);

// Leaf, root and wood litter lignin fractions
// Leaf and root fractions: Comins & McMurtrie 1993; Friend et al 1997
const double LIGCFRAC_LEAF = 0.2;
const double LIGCFRAC_ROOT = 0.16;
const double LIGCFRAC_WOOD = 0.3;	// TODO Check wood fraction

#endif // LPJ_GUESS_SOMDYNAM_H
