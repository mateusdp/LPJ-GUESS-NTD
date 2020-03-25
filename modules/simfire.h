///////////////////////////////////////////////////////////////////////////////////////
/// \file simfire.h
/// \brief SIMFIRE - SIMple FIRE module to compute burnt area  
///
/// \author Lars Nieradzik
/// $Date: 2015-08-25 09:19:28 +0200 (Tue, 25 Aug 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SIMFIRE_H
#define LPJ_GUESS_SIMFIRE_H

#include "guess.h"

// SIMFIRE biometypes
enum {SF_NOVEG, SF_CROP, SF_NEEDLELEAF, SF_BROADLEAF, SF_MIXED_FOREST, SF_SHRUBS, SF_SAVANNA, SF_TUNDRA, SF_BARREN};

/// Read SIMFIRE related data for a gridcell at beginning of simulation
/** Reads SIMFIRE relevant info from simfire_nput.bin:
 *  Hyde 3.1 population density
 *  Monthly fire climatology
 */
void getsimfiredata(Gridcell& gridcell);

/// Daily book-keeping for SIMFIRE-relevant variables
/** Updates SIMFIRE's Max Annual Mesterov Index
 *  and running mean of max annual FPAR (from canexch.cpp).
 *  Updates fire biome at beginning of the year.
 */
void simfire_accounting_gridcell(Gridcell& gridcell);

/// Calculate burned area in ha following Knorr 2014.
double simfire_burned_area(Gridcell& gridcell); 

#endif // LPJ_GUESS_SIMFIRE_H
