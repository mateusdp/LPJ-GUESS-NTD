///////////////////////////////////////////////////////////////////////////////////////
/// \file euroflux.cpp
/// \brief Extra code used by the Euroflux benchmarks
///
/// \author Joe Siltberg
/// $Date: 2015-11-13 16:25:45 +0100 (Fr, 13. Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "euroflux.h"
#include "parameters.h"
#include "guess.h"
#include "driver.h"

REGISTER_INPUT_MODULE("euroflux", EurofluxInput)
REGISTER_OUTPUT_MODULE("euroflux", EurofluxOutput)

using namespace GuessOutput;

// guess2008 - euroflux - new int to keep track of the simulation year
// Needed for management etc., used in vegetation dynamics
// century_year = 0, when date.year < nyear, i.e. during spin up. 
// century_year = 1, when date.year = nyear, i.e. 1901, 
// century_year = 80, when date.year = nyear+79, i.e. 1980, etc.
int century_year;

EurofluxData* current_stand_fluxdata = 0;

void EurofluxInput::init() {
	// First let base class initialize
	CRUInput::init();

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	bool eof = false;

	while (!eof) {

		EurofluxData edata;

		double dlon, dlat;
		xtring desc2;

		// New, local versions of these arrays
		xtring dom_spec[5];
		int dom_spec_dens[5];

		eof=!readfor(in_grid,"f,f,a;a;a;f;f;f;f;11i;f;i;i;5a;i;i;i;i;i",&dlon,&dlat,&edata.desc,
			&desc2,&edata.ver,&edata.tm,&edata.tc,&edata.pm,&edata.pc,
			edata.isfluxdata,&edata.soildepth,&edata.plantation_year,&edata.num_dominant_species,dom_spec,
			&dom_spec_dens[0],&dom_spec_dens[1],&dom_spec_dens[2],&dom_spec_dens[3],&dom_spec_dens[4]);


		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)

			int dsp;

			// Now read num_dominant_species lines from the gridlist files
			for (dsp = 0; dsp < edata.num_dominant_species; dsp++) {
				edata.dom_species[dsp] = dom_spec[dsp];
				edata.dom_species_density[dsp] = dom_spec_dens[dsp];
			}


			for (dsp = edata.num_dominant_species; dsp < 5; dsp++) {
				edata.dom_species[dsp] = "NONE";
				edata.dom_species_density[dsp] = 0;
			}

			eurofluxdata[std::make_pair(dlon, dlat)] = edata;
		}
	}

	fclose(in_grid);
}

namespace {

/// Adjusts some of the Soiltype objects members according to a site's soil depth
void adjust_soildepth(Soiltype& soiltype, double soildepth) {

	if (soildepth >= 600.0) {
		soiltype.awc[1] = min((soildepth-SOILDEPTH_UPPER),SOILDEPTH_LOWER) * (soiltype.awc[1]/SOILDEPTH_LOWER);
		soiltype.wp[1] = min((soildepth-SOILDEPTH_UPPER),SOILDEPTH_LOWER) * (soiltype.wp[1]/SOILDEPTH_LOWER);
		soiltype.wsats[1] = min((soildepth-SOILDEPTH_UPPER),SOILDEPTH_LOWER) * (soiltype.wsats[1]/SOILDEPTH_LOWER);
	}
	else {
		soiltype.awc[0] /= 2.0; // 25cm
		soiltype.awc[1] = min((soildepth-SOILDEPTH_UPPER/2.0),SOILDEPTH_LOWER) * (soiltype.awc[1]/SOILDEPTH_LOWER);
		soiltype.wp[0] /= 2.0;
		soiltype.wp[1] = min((soildepth-SOILDEPTH_UPPER/2.0),SOILDEPTH_LOWER) * (soiltype.wp[1]/SOILDEPTH_LOWER);
		soiltype.wsats[0] /= 2.0;
		soiltype.wsats[1] = min((soildepth-SOILDEPTH_UPPER/2.0),SOILDEPTH_LOWER) * (soiltype.wsats[1]/SOILDEPTH_LOWER);
	}

	soiltype.wtot = (soiltype.wtot/(SOILDEPTH_UPPER + SOILDEPTH_LOWER))*soildepth;
}

}

bool EurofluxInput::getgridcell(Gridcell& gridcell) {
	if (CRUInput::getgridcell(gridcell)) {
		
		// now give the soil depth too
		EurofluxData& edata = eurofluxdata[std::make_pair(gridcell.get_lon(), gridcell.get_lat())];
		adjust_soildepth(gridcell.soiltype, edata.soildepth);
		
		return true;
	}
	else {
		return false;
	}
}

bool EurofluxInput::getclimate(Gridcell& gridcell) {
	if (date.day == 0) {
		if (date.year < nyear_spinup) {
			century_year = 0;
		}
		else {
			century_year++;
		}
	}

	return CRUInput::getclimate(gridcell);
}

void EurofluxInput::adjust_raw_forcing_data(double lon,
                                            double lat,
                                            double hist_mtemp[NYEAR_HIST][12],
                                            double hist_mprec[NYEAR_HIST][12],
                                            double hist_msun[NYEAR_HIST][12]) {

	current_stand_fluxdata = &eurofluxdata[std::make_pair(lon, lat)];

	EurofluxData& coord = eurofluxdata[std::make_pair(lon, lat)];


	// *** Adjust all CRU temp and precip data to site conditions

	// Regression coefficients for this flux site, as read from the gridlist file
	double tempm = coord.tm;
	double tempc = coord.tc;
	double precipm = coord.pm;
	double precipc = coord.pc;

	for (int y = 0; y < NYEAR_HIST; y++) {
		for (int m = 0; m < 12; m++) {

			// Adjust CRU data to site conditions
			hist_mtemp[y][m] = tempm * hist_mtemp[y][m] + tempc;
			hist_mprec[y][m] = precipm * hist_mprec[y][m] + precipc;
				
			// Hack! Because negligible precipitation causes problems in the 
			// prdaily function (infinite loops). 
			if (hist_mprec[y][m] <= 1.0) hist_mprec[y][m] = 0.0;

		}
	}	

	// Get actual temp and precip data for the site, as well as NEE and latent heat flux

	// Create some strings 
	xtring fluxdirectory=param["flux_dir"].str;
	xtring fluxfilestart = "CEIP_EC_L4_m_";

	// Could possible get rid of the ver string, and try to open both v1 and v2...
	xtring fluxfileend = coord.ver;
	fluxfileend += ".txt";

	fluxdirectory += fluxfilestart;
	fluxdirectory+=coord.desc;
	fluxdirectory+="_";


	// Backup data in case NEE_st values are all -9999.0
	double NEE_or[NFLUXYEARS][12];
	double GPP_or[NFLUXYEARS][12];

	for (int y = 0; y < NFLUXYEARS; y++) {
		for (int m = 0; m < 12; m++) {
			NEE_or[y][m] = MISSING_DATA;
			GPP_or[y][m] = MISSING_DATA;
		}
	}


	// Extend array if we go beyond 2002.
	xtring fluxyears[NFLUXYEARS] = {"1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006"};

	int fyear = 0;

	// Loop from 1996 to 2006
	for (int y = NYEAR_HIST-NFLUXYEARS; y < NYEAR_HIST; y++) {

		// Is there flux data for this year? 
		if (coord.isfluxdata[fyear] == 1) {
				
			// Determine the full file name for this site and year
			xtring datafile = fluxdirectory;
			datafile += fluxyears[fyear];
			datafile += "_";
			datafile += fluxfileend;

			// test
			//datafile = fluxdirectory + "testin.txt";

			FILE* in_flux=fopen(datafile,"r");
			if (!in_flux) fail("getgridcell: could not open %s for input",(char*)datafile);

			bool eof = false;
			xtring header;

			// Read the header first. We don't use this.
			eof=!readfor(in_flux,"a",&header);
				
			int month = 0;

			// Latent heat of vapourisation [J/kg]
			const double LATENT_HEAT_VAP = 2510400.0; 

			// Seconds in a day
			const double SECS_IN_DAY = 24.0 * 60.0 * 60.0; 

			// Minimum quality required
			const double MIN_OBS_FREQ = 0.5;

			while (!eof) {
		
				bool useOriginalDataThisMonth = true; //_or data of filled data?

				// Each file has 13 rows and 30 columns.
				double sitedata[30];
				eof=!readfor(in_flux,"f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f",
				             &sitedata[0],&sitedata[1],&sitedata[2],&sitedata[3],&sitedata[4],&sitedata[5],&sitedata[6],
				             &sitedata[7],&sitedata[8],&sitedata[9],&sitedata[10],&sitedata[11],&sitedata[12],&sitedata[13],
				             &sitedata[14],&sitedata[15],&sitedata[16],&sitedata[17],&sitedata[18],&sitedata[19],&sitedata[20],
				             &sitedata[21],&sitedata[22],&sitedata[23],&sitedata[24],&sitedata[25],&sitedata[26],&sitedata[27],
				             &sitedata[28],&sitedata[29]);

				if (!eof) {

					// Read the relevant climate data from the site
					double mth			= sitedata[0];	// Month (1-12)
					double n_days		= sitedata[1];	// #days
					double Ta_f			= sitedata[4];	// degC
					double Ta_sqc		= sitedata[5];	// [0,1]
					double precip		= sitedata[10]; // mm day-1

					// Override the CRU data with the actual site climate data, where available.
					if (Ta_f != MISSING_DATA && Ta_sqc >= MIN_OBS_FREQ) // Only data of sufficient quality is stored
						hist_mtemp[y][month] = Ta_f;

					if (precip != MISSING_DATA)
						hist_mprec[y][month] = n_days*precip;
						

					// Now read the relevant soil and flux data from the site
					double swc			= sitedata[11]; // %vol
					double LE_f			= sitedata[14];	// W m-2 day-1
					double LE_sqc		= sitedata[15]; // [0,1]
					double NEE_st_fMDS	= sitedata[18];	// gC m-2 day-1
					double NEE_st_fMDSsqc = sitedata[19];	// [0,1]
					double GPP_st_MDS	= sitedata[20];	// gC m-2 day-1
					double NEE_or_fMDS	= sitedata[21];	// gC m-2 day-1
					double NEE_or_fMDSsqc = sitedata[22];	// [0,1]
					double GPP_or_MDS	= sitedata[23];	// gC m-2 day-1


					// Only AET data of sufficient quality is stored.
					// Convert from W m-2 day-1 to mm month-1
					if (LE_sqc >= MIN_OBS_FREQ)	
						current_stand_fluxdata->fluxAET[fyear][month] = n_days * SECS_IN_DAY / LATENT_HEAT_VAP * LE_f;

					// Only NEE_st data of sufficient quality is stored
					if (NEE_st_fMDSsqc >= MIN_OBS_FREQ && NEE_st_fMDS != MISSING_DATA) {
						useOriginalDataThisMonth = false; // No need replace this data with _or data	 
						current_stand_fluxdata->fluxNEE[fyear][month] = n_days * NEE_st_fMDS;
						current_stand_fluxdata->fluxGPP[fyear][month] = n_days * GPP_st_MDS;
					}

					// Back-up NEE_or data of sufficient quality
					if (NEE_or_fMDSsqc >= MIN_OBS_FREQ && NEE_or_fMDS != MISSING_DATA) {							
						NEE_or[fyear][month] = n_days * NEE_or_fMDS;
						GPP_or[fyear][month] = n_days * GPP_or_MDS;
					} else {
						NEE_or[fyear][month] = MISSING_DATA;
						GPP_or[fyear][month] = MISSING_DATA;
					}

						
					// Replace bad data with original data?
					if (useOriginalDataThisMonth) {
						current_stand_fluxdata->fluxNEE[fyear][month] = NEE_or[fyear][month];
						current_stand_fluxdata->fluxGPP[fyear][month] = GPP_or[fyear][month];
					}
				

					current_stand_fluxdata->fluxSWC[fyear][month] = swc;
				

					month++;

				} // if (!eof)

			} // while (!eof)


			fclose(in_flux);

			// Error?
			if (month != 12) {
				fail("\nError: could not read the data from the following flux file:\n%s\n", 
				     (char*)datafile);
			}

		} // isfluxdata
			
		fyear++;

	} // for
}



EurofluxOutput::EurofluxOutput() {
	// Files for EUROFLUX output
	declare_parameter("file_eurofluxmonthly_nee", &file_eurofluxmonthly_nee, 300, "EUROFLUX NEE monthly output file");
	declare_parameter("file_eurofluxmonthly_aet", &file_eurofluxmonthly_aet, 300, "EUROFLUX AET monthly output file");
	declare_parameter("file_eurofluxmonthly_gpp", &file_eurofluxmonthly_gpp, 300, "EUROFLUX GPP monthly output file");
	declare_parameter("file_eurofluxannual", &file_eurofluxannual, 300, "EUROFLUX annual output file");
	declare_parameter("file_eurofluxstats_nee", &file_eurofluxstats_nee, 300, "EUROFLUX NEE Statistics output file");
	declare_parameter("file_eurofluxstats_aet", &file_eurofluxstats_aet, 300, "EUROFLUX AET Statistics output file");
	declare_parameter("file_eurofluxstats_gpp", &file_eurofluxstats_gpp, 300, "EUROFLUX GPP Statistics output file");
}

void EurofluxOutput::init() {

	// create the output tables

	ColumnDescriptors month_columns_wide;
	xtring months[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
	for (int i = 0; i < 12; i++) {
		month_columns_wide += ColumnDescriptor(months[i], 10, 3);
	}

	ColumnDescriptors eurofluxannual_columns;
	xtring eurofluxannual_str[] = { "NEE_mod","NEE_obs","AET_mod","AET_obs","GPP_mod","GPP_obs",         // annual
	                                "sNEE_mod","sNEE_obs","sAET_mod","sAET_obs","sGPP_mod","sGPP_obs" }; // summer
	for (int i = 0; i < 12; i++) {
		eurofluxannual_columns += ColumnDescriptor(eurofluxannual_str[i], 10, 2);
	}

	ColumnDescriptors eurofluxstats_columns;
	xtring eurofluxstats_str[] = { "<OBS>","<MOD>","M","R","R2","EF","CD","RMSE","RMSE_std","t" };

	eurofluxstats_columns += ColumnDescriptor("N", 10, 0);
	for (int i = 0; i < 10; i++) {
		eurofluxstats_columns += ColumnDescriptor(eurofluxstats_str[i], 10, 2);
	}

	create_output_table(out_eurofluxannual, file_eurofluxannual, eurofluxannual_columns);
	create_output_table(out_eurofluxmonthly_nee, file_eurofluxmonthly_nee, month_columns_wide);
	create_output_table(out_eurofluxmonthly_aet, file_eurofluxmonthly_aet, month_columns_wide);
	create_output_table(out_eurofluxmonthly_gpp, file_eurofluxmonthly_gpp, month_columns_wide);
	create_output_table(out_eurofluxstats_nee, file_eurofluxstats_nee, eurofluxstats_columns);
	create_output_table(out_eurofluxstats_aet, file_eurofluxstats_aet, eurofluxstats_columns);
	create_output_table(out_eurofluxstats_gpp, file_eurofluxstats_gpp, eurofluxstats_columns);
}

///////////////////////////////////////////////////////////////////////////////////////
// calculateAnnualFluxSums - euroflux
// Called by outannual at the end of each flux year
void calculateAnnualFluxSums(const int yr, 
								double& annNEE_obs, double& annNEE_mod, double &sumNEE_obs, double& sumNEE_mod, 
								double& annAET_obs, double& annAET_mod, double &sumAET_obs, double& sumAET_mod, 
								double& annGPP_obs, double& annGPP_mod, double &sumGPP_obs, double& sumGPP_mod) {

	int mth;

	// Initialise to 0.0;

	annNEE_obs = 0.0;
	annNEE_mod = 0.0;
	sumNEE_obs = 0.0;
	sumNEE_mod = 0.0;
	
	annAET_obs = 0.0;
	annAET_mod = 0.0;
	sumAET_obs = 0.0;
	sumAET_mod = 0.0;

	annGPP_obs = 0.0;
	annGPP_mod = 0.0;
	sumGPP_obs = 0.0;
	sumGPP_mod = 0.0;


	int yrNEE_obs = 0;
	int yrAET_obs = 0;
	int yrGPP_obs = 0;

	int jjaNEE_obs = 0;
	int jjaAET_obs = 0;
	int jjaGPP_obs = 0;

	for (mth = 0; mth < 12; mth++) {

		if (current_stand_fluxdata->fluxNEE[yr][mth] != MISSING_DATA) {
			
			// Annual stats
			yrNEE_obs++;
			annNEE_obs += current_stand_fluxdata->fluxNEE[yr][mth];
			annNEE_mod += current_stand_fluxdata->modelNEE[yr][mth];
			
			// Summer (JJA) stats
			if (mth >= 6 && mth <= 8) {
				jjaNEE_obs++;
				sumNEE_obs += current_stand_fluxdata->fluxNEE[yr][mth];
				sumNEE_mod += current_stand_fluxdata->modelNEE[yr][mth];
			}

		} // NEE

		if (current_stand_fluxdata->fluxAET[yr][mth] != MISSING_DATA) {

			// Annual stats
			yrAET_obs++;
			annAET_obs += current_stand_fluxdata->fluxAET[yr][mth];
			annAET_mod += current_stand_fluxdata->modelAET[yr][mth];
			
			// Summer (JJA) stats
			if (mth >= 6 && mth <= 8) {
				jjaAET_obs++;
				sumAET_obs += current_stand_fluxdata->fluxAET[yr][mth];
				sumAET_mod += current_stand_fluxdata->modelAET[yr][mth];
			}

		} // AET

		if (current_stand_fluxdata->fluxGPP[yr][mth] != MISSING_DATA) {
		
			// Annual stats
			yrGPP_obs++;
			annGPP_obs += current_stand_fluxdata->fluxGPP[yr][mth];
			annGPP_mod += current_stand_fluxdata->modelGPP[yr][mth];
			
			// Summer (JJA) stats
			if (mth >= 6 && mth <= 8) {
				jjaGPP_obs++;
				sumGPP_obs += current_stand_fluxdata->fluxGPP[yr][mth];
				sumGPP_mod += current_stand_fluxdata->modelGPP[yr][mth];
			}
		
		} // GPP

	} // mth loop


	// Data?
	if (yrNEE_obs == 0) {
		annNEE_obs = MISSING_DATA;
		annNEE_mod = MISSING_DATA;
	}

	if (jjaNEE_obs == 0) {
		sumNEE_obs = MISSING_DATA;
		sumNEE_mod = MISSING_DATA;
	}

	if (yrAET_obs == 0) {
		annAET_obs = MISSING_DATA;
		annAET_mod = MISSING_DATA;
	}

	if (jjaAET_obs == 0) {
		sumAET_obs = MISSING_DATA;
		sumAET_mod = MISSING_DATA;
	}

	if (yrGPP_obs == 0) {
		annGPP_obs = MISSING_DATA;
		annGPP_mod = MISSING_DATA;
	}

	if (jjaGPP_obs == 0) {
		sumGPP_obs = MISSING_DATA;
		sumGPP_mod = MISSING_DATA;
	}


}



///////////////////////////////////////////////////////////////////////////////////////
// calculateEurofluxStats - euroflux
// Called by outannual at the end of the last day of the last simulation year

void calculateEurofluxStats(double lon,
                            double lat,
                            Table out_stats_nee, 
                            Table out_stats_aet, 
                            Table out_stats_gpp) {

	int yr, mth;

	// First calculate the modelled and observed averages for NEE, AET and GPP
	// ***********************************************************************

	double meanNEE_obs = 0.0;
	double meanNEE_mod = 0.0;
	int numNEE_obs = 0;

	double meanAET_obs = 0.0;
	double meanAET_mod = 0.0;
	int numAET_obs = 0;
	
	double meanGPP_obs = 0.0;
	double meanGPP_mod = 0.0;
	int numGPP_obs = 0;

	
	for (yr = 0; yr < NFLUXYEARS; yr++) {

		for (int mth = 0; mth < 12; mth++) {

			if (current_stand_fluxdata->fluxNEE[yr][mth] != MISSING_DATA) {
				
				// Overall stats
				numNEE_obs++;
				meanNEE_obs += current_stand_fluxdata->fluxNEE[yr][mth];
				meanNEE_mod += current_stand_fluxdata->modelNEE[yr][mth];
	
			} // NEE

			if (current_stand_fluxdata->fluxAET[yr][mth] != MISSING_DATA) {
				numAET_obs++;
				meanAET_obs += current_stand_fluxdata->fluxAET[yr][mth];
				meanAET_mod += current_stand_fluxdata->modelAET[yr][mth];

			} // AET

			if (current_stand_fluxdata->fluxGPP[yr][mth] != MISSING_DATA) {
				numGPP_obs++;
				meanGPP_obs += current_stand_fluxdata->fluxGPP[yr][mth];
				meanGPP_mod += current_stand_fluxdata->modelGPP[yr][mth];
						
			} // GPP

		} // mth loop


	} // yr loop


	// Overall averages
	meanNEE_obs /= numNEE_obs;
	meanAET_obs /= numAET_obs;
	meanGPP_obs /= numGPP_obs;

	meanNEE_mod /= numNEE_obs;
	meanAET_mod /= numAET_obs;
	meanGPP_mod /= numGPP_obs;




	// Now prepare to calculate the correlation coefficients, RMSE, M, and Sd 
	// statistics for NEE, AET and GPP
	// ***********************************************************************


	double r_NEE = 0.0;
	double rmse_NEE = 0.0;
	double rmse_NEE_std = 0.0;
	double m_NEE = 0.0;
	double sd_NEE_sq = 0.0;
	double t_NEE = 0.0;

	double r_AET = 0.0;
	double rmse_AET = 0.0;
	double rmse_AET_std = 0.0;
	double m_AET = 0.0;
	double sd_AET_sq = 0.0;
	double t_AET= 0.0;

	double r_GPP = 0.0;
	double rmse_GPP = 0.0;
	double rmse_GPP_std = 0.0;
	double m_GPP = 0.0;
	double sd_GPP_sq = 0.0;
	double t_GPP = 0.0;

	double NEE_obs_diff = 0.0; // SUM(O-<O>)
	double NEE_mod_diff = 0.0; // SUM(P-<P>)
	double NEE_obs_sqdiff = 0.0; // SUM((O-<O>)(O-<O>))
	double NEE_mod_sqdiff = 0.0; // SUM((P-<P>)(P-<P>))
	double NEE_obs_minus_mod = 0.0; // SUM(O-P)
	double NEE_obs_minus_mod_sq = 0.0; // SUM((O-P)(O-P))
	double NEE_obs_diff_times_mod_diff = 0.0; // SUM((O-<O>)(P-<P>))

	double AET_obs_diff = 0.0; // SUM(O-<O>)
	double AET_mod_diff = 0.0; // SUM(P-<P>)
	double AET_obs_sqdiff = 0.0; // SUM((O-<O>)(O-<O>))
	double AET_mod_sqdiff = 0.0; // SUM((P-<P>)(P-<P>))
	double AET_obs_minus_mod = 0.0; // SUM(O-P)
	double AET_obs_minus_mod_sq = 0.0; // SUM((O-P)(O-P))
	double AET_obs_diff_times_mod_diff = 0.0; // SUM((O-<O>)(P-<P>))

	double GPP_obs_diff = 0.0; // SUM(O-<O>)
	double GPP_mod_diff = 0.0; // SUM(P-<P>)
	double GPP_obs_sqdiff = 0.0; // SUM((O-<O>)(O-<O>))
	double GPP_mod_sqdiff = 0.0; // SUM((P-<P>)(P-<P>))
	double GPP_obs_minus_mod = 0.0; // SUM(O-P)
	double GPP_obs_minus_mod_sq = 0.0; // SUM((O-P)(O-P))
	double GPP_obs_diff_times_mod_diff = 0.0; // SUM((O-<O>)(P-<P>))


	for (yr = 0; yr < NFLUXYEARS; yr++) {
		for (int mth = 0; mth < 12; mth++) {

			// NEE stats
			if (current_stand_fluxdata->fluxNEE[yr][mth] != MISSING_DATA) {
				
				NEE_obs_diff += current_stand_fluxdata->fluxNEE[yr][mth] - meanNEE_obs;
				NEE_obs_sqdiff += (current_stand_fluxdata->fluxNEE[yr][mth] - meanNEE_obs) * (current_stand_fluxdata->fluxNEE[yr][mth] - meanNEE_obs);

				NEE_mod_diff += current_stand_fluxdata->modelNEE[yr][mth] - meanNEE_mod;
				NEE_mod_sqdiff += (current_stand_fluxdata->modelNEE[yr][mth] - meanNEE_mod) * (current_stand_fluxdata->modelNEE[yr][mth] - meanNEE_mod);

				NEE_obs_minus_mod += current_stand_fluxdata->fluxNEE[yr][mth] - current_stand_fluxdata->modelNEE[yr][mth];
				NEE_obs_minus_mod_sq += (current_stand_fluxdata->fluxNEE[yr][mth] - current_stand_fluxdata->modelNEE[yr][mth]) * (current_stand_fluxdata->fluxNEE[yr][mth] - current_stand_fluxdata->modelNEE[yr][mth]);
				
				NEE_obs_diff_times_mod_diff += (current_stand_fluxdata->fluxNEE[yr][mth] - meanNEE_obs) * (current_stand_fluxdata->modelNEE[yr][mth] - meanNEE_mod);

			}

			// AET stats
			if (current_stand_fluxdata->fluxAET[yr][mth] != MISSING_DATA) {

				AET_obs_diff += current_stand_fluxdata->fluxAET[yr][mth] - meanAET_obs;
				AET_obs_sqdiff += (current_stand_fluxdata->fluxAET[yr][mth] - meanAET_obs) * (current_stand_fluxdata->fluxAET[yr][mth] - meanAET_obs);

				AET_mod_diff += current_stand_fluxdata->modelAET[yr][mth] - meanAET_mod;
				AET_mod_sqdiff += (current_stand_fluxdata->modelAET[yr][mth] - meanAET_mod) * (current_stand_fluxdata->modelAET[yr][mth] - meanAET_mod);

				AET_obs_minus_mod += current_stand_fluxdata->fluxAET[yr][mth] - current_stand_fluxdata->modelAET[yr][mth];
				AET_obs_minus_mod_sq += (current_stand_fluxdata->fluxAET[yr][mth] - current_stand_fluxdata->modelAET[yr][mth]) * (current_stand_fluxdata->fluxAET[yr][mth] - current_stand_fluxdata->modelAET[yr][mth]);

				AET_obs_diff_times_mod_diff += (current_stand_fluxdata->fluxAET[yr][mth] - meanAET_obs) * (current_stand_fluxdata->modelAET[yr][mth] - meanAET_mod);

			}
			
			// GPP stats
			if (current_stand_fluxdata->fluxGPP[yr][mth] != MISSING_DATA) {

				GPP_obs_diff += current_stand_fluxdata->fluxGPP[yr][mth] - meanGPP_obs;
				GPP_obs_sqdiff += (current_stand_fluxdata->fluxGPP[yr][mth] - meanGPP_obs) * (current_stand_fluxdata->fluxGPP[yr][mth] - meanGPP_obs);

				GPP_mod_diff += current_stand_fluxdata->modelGPP[yr][mth] - meanGPP_mod;
				GPP_mod_sqdiff += (current_stand_fluxdata->modelGPP[yr][mth] - meanGPP_mod) * (current_stand_fluxdata->modelGPP[yr][mth] - meanGPP_mod);

				GPP_obs_minus_mod += current_stand_fluxdata->fluxGPP[yr][mth] - current_stand_fluxdata->modelGPP[yr][mth];
				GPP_obs_minus_mod_sq += (current_stand_fluxdata->fluxGPP[yr][mth] - current_stand_fluxdata->modelGPP[yr][mth]) * (current_stand_fluxdata->fluxGPP[yr][mth] - current_stand_fluxdata->modelGPP[yr][mth]);

				GPP_obs_diff_times_mod_diff += (current_stand_fluxdata->fluxGPP[yr][mth] - meanGPP_obs) * (current_stand_fluxdata->modelGPP[yr][mth] - meanGPP_mod);

			}
		} // for mth
	} // for yr


	// *** NEE stats ***
	// ***********************************************************************

	rmse_NEE = (100.0 / meanNEE_obs) * sqrt(NEE_obs_minus_mod_sq / numNEE_obs);
	rmse_NEE_std = sqrt(NEE_obs_minus_mod_sq / numNEE_obs);
	m_NEE = NEE_obs_minus_mod / numNEE_obs;
	r_NEE = NEE_obs_diff_times_mod_diff / sqrt(NEE_obs_sqdiff * NEE_mod_sqdiff);

	for (yr = 0; yr < NFLUXYEARS; yr++) {
		for (int mth = 0; mth < 12; mth++) {

			// sd_NEE_sq
			if (current_stand_fluxdata->fluxNEE[yr][mth] != MISSING_DATA) {		
				sd_NEE_sq += (current_stand_fluxdata->fluxNEE[yr][mth] - current_stand_fluxdata->modelNEE[yr][mth] - m_NEE) *
				(current_stand_fluxdata->fluxNEE[yr][mth] - current_stand_fluxdata->modelNEE[yr][mth] - m_NEE);
			}

		} // for mth
	} // for yr

	sd_NEE_sq /= (numNEE_obs-1);
	t_NEE = m_NEE * sqrt(numNEE_obs) / sqrt(sd_NEE_sq);
	

	// *** AET stats ***
	// ***********************************************************************

	rmse_AET = (100.0 / meanAET_obs) * sqrt(AET_obs_minus_mod_sq / numAET_obs);
	rmse_AET_std = sqrt(AET_obs_minus_mod_sq / numAET_obs);
	m_AET = AET_obs_minus_mod / numAET_obs;
	r_AET = AET_obs_diff_times_mod_diff / sqrt(AET_obs_sqdiff * AET_mod_sqdiff);

	for (yr = 0; yr < NFLUXYEARS; yr++) {
		for (int mth = 0; mth < 12; mth++) {

			// sd_AET_sq
			if (current_stand_fluxdata->fluxAET[yr][mth] != MISSING_DATA) {		
				sd_AET_sq += (current_stand_fluxdata->fluxAET[yr][mth] - current_stand_fluxdata->modelAET[yr][mth] - m_AET) *
				(current_stand_fluxdata->fluxAET[yr][mth] - current_stand_fluxdata->modelAET[yr][mth] - m_AET);
			}

		} // for mth
	} // for yr

	sd_AET_sq /= (numAET_obs-1);
	t_AET = m_AET * sqrt(numAET_obs) / sqrt(sd_AET_sq); 
	

	// *** GPP stats ***
	// ***********************************************************************

	rmse_GPP = (100.0 / meanGPP_obs) * sqrt(GPP_obs_minus_mod_sq / numGPP_obs);
	rmse_GPP_std = sqrt(GPP_obs_minus_mod_sq / numGPP_obs);
	m_GPP = GPP_obs_minus_mod / numGPP_obs;
	r_GPP = GPP_obs_diff_times_mod_diff / sqrt(GPP_obs_sqdiff * GPP_mod_sqdiff);

	for (yr = 0; yr < NFLUXYEARS; yr++) {
		for (int mth = 0; mth < 12; mth++) {

			// sd_GPP_sq
			if (current_stand_fluxdata->fluxGPP[yr][mth] != MISSING_DATA) {		
				sd_GPP_sq += (current_stand_fluxdata->fluxGPP[yr][mth] - current_stand_fluxdata->modelGPP[yr][mth] - m_GPP) *
				(current_stand_fluxdata->fluxGPP[yr][mth] - current_stand_fluxdata->modelGPP[yr][mth] - m_GPP);
			}

		} // for mth
	} // for yr

	sd_GPP_sq /= (numGPP_obs-1);
	t_GPP = m_GPP * sqrt(numGPP_obs) / sqrt(sd_GPP_sq); 


	// OUTPUT stats
	// ***********************************************************************

	OutputRows out(output_channel, lon, lat, date.get_calendar_year());

	double ef = (NEE_obs_sqdiff - NEE_obs_minus_mod_sq)/NEE_obs_sqdiff; // Modelling Efficiency
	double cd = NEE_obs_sqdiff/NEE_obs_minus_mod_sq; // Coefficient of determination

	out.add_value(out_stats_nee, numNEE_obs);
	out.add_value(out_stats_nee, meanNEE_obs);
	out.add_value(out_stats_nee, meanNEE_mod);
	out.add_value(out_stats_nee, m_NEE);
	out.add_value(out_stats_nee, r_NEE);
	out.add_value(out_stats_nee, r_NEE*r_NEE);
	out.add_value(out_stats_nee, ef);
	out.add_value(out_stats_nee, cd);
	out.add_value(out_stats_nee, rmse_NEE);
	out.add_value(out_stats_nee, rmse_NEE_std);
	out.add_value(out_stats_nee, t_NEE);
	
	ef = (AET_obs_sqdiff - AET_obs_minus_mod_sq)/AET_obs_sqdiff; // Modelling Efficiency
	cd = AET_obs_sqdiff/AET_obs_minus_mod_sq; // Coefficient of determination

	out.add_value(out_stats_aet, numAET_obs);
	out.add_value(out_stats_aet, meanAET_obs);
	out.add_value(out_stats_aet, meanAET_mod);
	out.add_value(out_stats_aet, m_AET);
	out.add_value(out_stats_aet, r_AET);
	out.add_value(out_stats_aet, r_AET*r_AET);
	out.add_value(out_stats_aet, ef);
	out.add_value(out_stats_aet, cd);
	out.add_value(out_stats_aet, rmse_AET);
	out.add_value(out_stats_aet, rmse_AET_std);
	out.add_value(out_stats_aet, t_AET);

	ef = (GPP_obs_sqdiff - GPP_obs_minus_mod_sq)/GPP_obs_sqdiff; // Modelling Efficiency
	cd = GPP_obs_sqdiff/GPP_obs_minus_mod_sq; // Coefficient of determination

	out.add_value(out_stats_gpp, numGPP_obs);
	out.add_value(out_stats_gpp, meanGPP_obs);
	out.add_value(out_stats_gpp, meanGPP_mod);
	out.add_value(out_stats_gpp, m_GPP);
	out.add_value(out_stats_gpp, r_GPP);
	out.add_value(out_stats_gpp, r_GPP*r_GPP);
	out.add_value(out_stats_gpp, ef);
	out.add_value(out_stats_gpp, cd);
	out.add_value(out_stats_gpp, rmse_GPP);
	out.add_value(out_stats_gpp, rmse_GPP_std);
	out.add_value(out_stats_gpp, t_GPP);
}

void EurofluxOutput::outannual(Gridcell& gridcell) {

	// flux years are 1996 - 2002
	if (date.year >= nyear_spinup + CRUInput::NYEAR_HIST - NFLUXYEARS &&
	    date.year <= nyear_spinup + CRUInput::NYEAR_HIST - 1) {

		OutputRows out(output_channel, gridcell.get_lon(), gridcell.get_lat(), date.get_calendar_year());

		std::vector<double> mgpp(12);
		std::vector<double> maet(12);
		std::vector<double> mnee(12);

		Gridcell::iterator gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {

			Stand& stand = *gc_itr;
			stand.firstobj();
			
			while (stand.isobj) {
				Patch& patch = stand.getobj();

				double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();
			
				for (int m = 0; m < 12; ++m) {
					double gpp = patch.fluxes.get_monthly_flux(Fluxes::GPP, m);
					double ra = patch.fluxes.get_monthly_flux(Fluxes::RA, m);
					double rh = patch.fluxes.get_monthly_flux(Fluxes::SOILC, m);
					double npp = gpp - ra;
					mnee[m] += (rh - npp)*to_gridcell_average;
					mgpp[m] += gpp*to_gridcell_average;
					maet[m] += patch.maet[m]*to_gridcell_average;
				}
				
				stand.nextobj();
			}

			++gc_itr;
		}

		// euroflux
		// Print monthly EUROFLUX data
		for (int m = 0; m < 12; m++) {

			 int year = date.year-nyear_spinup-95;
			 double flux_nee = current_stand_fluxdata->fluxNEE[year][m];
			 double flux_aet = current_stand_fluxdata->fluxAET[year][m];
			 double flux_gpp = current_stand_fluxdata->fluxGPP[year][m];

			 // convert to the same units as our regular mnee, maet and mgpp files
			 if (flux_nee != MISSING_DATA) {
				  flux_nee *= 0.001;
			 }

			 if (flux_gpp != MISSING_DATA) {
				  flux_gpp *= 0.001;
			 }

			 out.add_value(out_eurofluxmonthly_nee, flux_nee);
			 out.add_value(out_eurofluxmonthly_aet, flux_aet);
			 out.add_value(out_eurofluxmonthly_gpp, flux_gpp);

			 // Save modelled flux data
			 current_stand_fluxdata->modelNEE[year][m] = -1000.0*mnee[m];  // gC/m2/month
			 current_stand_fluxdata->modelAET[year][m] = maet[m];			   // mm/month
			 current_stand_fluxdata->modelGPP[year][m] = 1000.0*mgpp[m];   // gC/m2/month
		}


		// Calculate and print yearly EUROFLUX data
			
		// Variables for annual and summer stats
		double annNEE_obs, annNEE_mod, sumNEE_obs, sumNEE_mod;
		double annAET_obs, annAET_mod, sumAET_obs, sumAET_mod;
		double annGPP_obs, annGPP_mod, sumGPP_obs, sumGPP_mod;

		calculateAnnualFluxSums(date.year-nyear_spinup-95, 
								annNEE_obs, annNEE_mod, sumNEE_obs, sumNEE_mod, 
								annAET_obs, annAET_mod, sumAET_obs, sumAET_mod, 
								annGPP_obs, annGPP_mod, sumGPP_obs, sumGPP_mod);

		// annual
		out.add_value(out_eurofluxannual, annNEE_mod);
		out.add_value(out_eurofluxannual, annNEE_obs);
		out.add_value(out_eurofluxannual, annAET_mod);
		out.add_value(out_eurofluxannual, annAET_obs);
		out.add_value(out_eurofluxannual, annGPP_mod);
		out.add_value(out_eurofluxannual, annGPP_obs);
		// summer
		out.add_value(out_eurofluxannual, sumNEE_mod);
		out.add_value(out_eurofluxannual, sumNEE_obs);
		out.add_value(out_eurofluxannual, sumAET_mod);
		out.add_value(out_eurofluxannual, sumAET_obs);
		out.add_value(out_eurofluxannual, sumGPP_mod);
		out.add_value(out_eurofluxannual, sumGPP_obs);
	}

	if (date.year == nyear_spinup + CRUInput::NYEAR_HIST - 1) {
		calculateEurofluxStats(gridcell.get_lon(),
		                       gridcell.get_lat(),
		                       out_eurofluxstats_nee,
		                       out_eurofluxstats_aet,
		                       out_eurofluxstats_gpp);
	}
	
}

void EurofluxOutput::outdaily(Gridcell& gridcell) {

}
