///////////////////////////////////////////////////////////////////////////////////////
/// \file growth.h
/// \brief The growth module header file
///
/// Vegetation C allocation, litter production, tissue turnover
/// leaf phenology, allometry and growth.
///
/// \author Ben Smith
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GROWTH_H
#define LPJ_GUESS_GROWTH_H

#include "guess.h"

double fracmass_lpj(double fpc_low,double fpc_high,Individual& indiv);
void leaf_phenology(Patch& patch,Climate& climate);
bool allometry(Individual& indiv); // guess2008 - now returns bool instead of void
void allocation_init(double bminit,double ltor,Individual& indiv);
void growth(Stand& stand,Patch& patch);
void growth_natural_daily(Stand& stand, Patch& patch);
void lai_update(Patch& patch);
//void turnover(double turnover_leaf, double turnover_root, double turnover_sap,
//	lifeformtype lifeform, landcovertype landcover, double& cmass_leaf, double& cmass_root, double& cmass_sap,
//	double& cmass_heart, double& nmass_leaf, double& nmass_root, double& nmass_sap,
//	double& nmass_heart, double& litter_leaf, double& litter_root,
//	double& nmass_litter_leaf, double& nmass_litter_root,
//	double& longterm_nstore, double &max_n_storage,
//	bool alive);
void turnover_np(double turnover_leaf, double turnover_root, double turnover_sap, double nreloc_ind, double preloc_ind,
	lifeformtype lifeform, landcovertype landcover, double& cmass_leaf, double& cmass_root, double& cmass_myco, double& cmass_sap,
	double& cmass_heart, double& nmass_leaf, double& nmass_root, double& nmass_sap,
	double& nmass_heart, double& pmass_leaf, double& pmass_root, double& pmass_sap,
	double& pmass_heart, double& litter_leaf, double& litter_root, double& litter_myco, double& cmass_leaf_root_turnover,
	double& nmass_litter_leaf, double& nmass_litter_root,
	double& pmass_litter_leaf, double& pmass_litter_root,
	double& longterm_nstore, double &max_n_storage,
	double& longterm_pstore, double &max_p_storage,
	bool alive);

#endif // LPJ_GUESS_GROWTH_H
