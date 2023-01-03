///////////////////////////////////////////////////////////////////////////////////////
/// \file driver.h
/// \brief Environmental driver calculation/transformation
///
/// \author Ben Smith
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_DRIVER_H
#define LPJ_GUESS_DRIVER_H

#include "guess.h"
#include <limits>
#include <map>
double randfrac(long& seed);

void interp_monthly_means_conserve(const double* mvals, double* dvals,
	double minimum = -std::numeric_limits<double>::max(),
	double maximum = std::numeric_limits<double>::max());
void interp_monthly_totals_conserve(const double* mvals, double* dvals,
	double minimum = -std::numeric_limits<double>::max(),
	double maximum = std::numeric_limits<double>::max());
void distribute_ndep(const double* mNH4dry, const double* mNO3dry, 
	const double* mNH4wet, const double* mNO3wet,
	const double* dprec, double* dNH4dep, double* dNO3dep);
void prdaily(double* mval_prec, double* dval_prec, double* mval_wet, long& seed, bool truncate = true);
void dailyaccounting_gridcell(Gridcell& gridcell);
void dailyaccounting_stand(Stand& stand);
void dailyaccounting_patch(Patch& patch);
void respiration_temperature_response(double temp,double& gtemp);
void daylengthinsoleet(Climate& climate);
double get_local_temp(Gridcell& gridcell, double gridcell_temp, double gridcell_elevation, double local_elevation);
void climate_seasonality(Gridcell& gridcell);

#endif // LPJ_GUESS_DRIVER_H
