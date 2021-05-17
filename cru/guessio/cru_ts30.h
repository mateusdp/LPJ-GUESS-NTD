///////////////////////////////////////////////////////////////////////////////////////
/// \file cru_ts30.h
/// \brief Functions for reading the CRU-NCEP data set from binary FastArchive format.
///
/// The binary files contain CRU-NCEP half-degree global historical climate data
/// for 1901-2015.
///
/// $Date: 2019-11-06 17:45:44 +0100 (Wed, 06 Nov 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CRU_TS30_H
#define LPJ_GUESS_CRU_TS30_H

namespace CRU_FastArchive {

/// number of years of historical climate
/** CRU-NCEP v7 has 115 years of data (1901-2015) */
const int NYEAR_HIST=115;

/// calendar year corresponding to first year in CRU climate data set
static const int FIRSTHISTYEAR=1901;

/// Determine temp, precip, sunshine & soilcode
bool searchcru(char* cruark,double dlon,double dlat,int& soilcode,
				double mtemp[NYEAR_HIST][12],
				double mprec[NYEAR_HIST][12],
				double msun[NYEAR_HIST][12]);

/// Determine elevation, frs frq, wet frq & DTR
bool searchcru_misc(char* cruark,double dlon,double dlat,int& elevation,
					double mfrs[NYEAR_HIST][12],
					double mwet[NYEAR_HIST][12],
					double mdtr[NYEAR_HIST][12],
					double mwind[NYEAR_HIST][12],
					double mrhum[NYEAR_HIST][12]);

/// Returns CRU data from the nearest cell to (lon,lat) within a given search radius
/** lon and lat are set to the coordinates of the found CRU gridcell, if found */
bool findnearestCRUdata(double searchradius, char* cruark, double& lon, double& lat, int& scode, 
						double hist_mtemp1[NYEAR_HIST][12], 
						double hist_mprec1[NYEAR_HIST][12], 
						double hist_msun1[NYEAR_HIST][12]);
}

#endif // LPJ_GUESS_CRU_TS30_H
