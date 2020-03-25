///////////////////////////////////////////////////////////////////////////////////////
/// \file weathergen.h
/// \brief Global Weather GENerator 
///
/// \author Lars Nieradzik
/// $Date: 2017-11-24 15:04:09 +0200 (Fri, 24 Nov 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "guess.h"

#ifndef WEATHERGEN_H
#define WEATHERGEN_H

/// GWGen - Global Weather GENerator
/** 
 * A weathergenerator for the use with e.g. BLAZE when wind and/or rel. 
 * humidity is needed. 
 */
void weathergen_get_met(Gridcell& gridcell, double* in_mtemp, double* in_mprec, 
		   double* in_mwetm, double* in_msol, double* in_mdtr, 
		   double* in_mwind, double* in_rhum, double* out_temp,
		   double* out_dprec,double* out_dsol,double* out_ddtr,
		   double* out_dwind,double* out_rhum);
	
#endif
