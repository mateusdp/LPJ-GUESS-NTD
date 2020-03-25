///////////////////////////////////////////////////////////////////////////////////////
/// \file euroflux.h
/// \brief Extra code used by the Euroflux benchmarks
///
/// \author Joe Siltberg
/// $Date: 2015-11-13 16:25:45 +0100 (Fr, 13. Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_EUROFLUX_H
#define LPJ_GUESS_EUROFLUX_H

#include "cruinput.h"
#include "outputmodule.h"
#include <gutil.h>

/// The value used for missing data in the EUROFLUX files.
const double MISSING_DATA = -9999.0;

/// The number of EUROFLUX years, currently 1996-2006, inclusive.
const int NFLUXYEARS=11;

/// Type for storing EUROFLUX grid cell information
struct EurofluxData {

	xtring desc; 

	xtring ver;
	double tm;
	double tc;
	double pm;
	double pc;

	int isfluxdata[NFLUXYEARS];

	double soildepth;
	int plantation_year;
	int num_dominant_species;
	xtring dom_species[5];
	int dom_species_density[5];
	
	// Flux data for a site
	double fluxNEE[NFLUXYEARS][12];
	double fluxAET[NFLUXYEARS][12];
	double fluxGPP[NFLUXYEARS][12];
	double fluxSWC[NFLUXYEARS][12];

	// Modelled flux data for the same site
	double modelNEE[NFLUXYEARS][12];
	double modelAET[NFLUXYEARS][12];
	double modelGPP[NFLUXYEARS][12];

	EurofluxData() {
		// initialise EUROFLUX arrays with missing values; 
		for (int yr = 0; yr < NFLUXYEARS; yr++) {
			for (int mth = 0; mth < 12; mth++) {
				fluxNEE[yr][mth] = MISSING_DATA;
				fluxAET[yr][mth] = MISSING_DATA;
				fluxGPP[yr][mth] = MISSING_DATA;
				fluxSWC[yr][mth] = MISSING_DATA;
			}
		}

		// New initialisation
		plantation_year = -1;
		num_dominant_species = 0;
		desc = ""; 

	}
};

/// Input module for EUROFLUX benchmark
/** This is a subclass of the CRU input module. The subclass
 *  will alter the CRU forcing data according to site data,
 *  and also update the soiltype with soildepth information
 *  after the base class has initialized it according to soil code.
 */
class EurofluxInput : public CRUInput {
public:
	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

protected:
	void adjust_raw_forcing_data(double lon,
	                             double lat,
	                             double hist_mtemp[NYEAR_HIST][12],
	                             double hist_mprec[NYEAR_HIST][12],
	                             double hist_msun[NYEAR_HIST][12]);

private:
	std::map<std::pair<double, double>, EurofluxData> eurofluxdata;
};


/// Output module for the extra files for EUROFLUX
class EurofluxOutput : public GuessOutput::OutputModule {
public:
	EurofluxOutput();

	void init();

	void outannual(Gridcell& gridcell);

	void outdaily(Gridcell& gridcell);
	void openlocalfiles(Gridcell& gridcell) {};
	void closelocalfiles(Gridcell& gridcell) {};

private:
	// Files for EUROFLUX output and stats
	xtring file_eurofluxannual;
	xtring file_eurofluxmonthly_nee, file_eurofluxmonthly_aet, file_eurofluxmonthly_gpp;
	xtring file_eurofluxstats_nee, file_eurofluxstats_aet, file_eurofluxstats_gpp;

	// Output tables
	GuessOutput::Table out_eurofluxannual;
	GuessOutput::Table out_eurofluxmonthly_nee, out_eurofluxmonthly_aet, out_eurofluxmonthly_gpp;
	GuessOutput::Table out_eurofluxstats_nee, out_eurofluxstats_aet, out_eurofluxstats_gpp;
};

#endif // LPJ_GUESS_EUROFLUX_H
