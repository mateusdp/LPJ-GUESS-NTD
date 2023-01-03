////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cropsowing.h
/// \brief Seasonality and sowing date calculations
/// \author Mats Lindeskog
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CROPSOWING_H
#define LPJ_GUESS_CROPSOWING_H

/// Monitors climate history relevant for sowing date calculation. Calculates initial sowing dates/windows
void crop_sowing_gridcell(Gridcell& gridcell);
/// Handles sowing date calculations for crop pft:s on patch level
void crop_sowing_patch(Patch& patch);
#endif // LPJ_GUESS_CROPSOWING_H
