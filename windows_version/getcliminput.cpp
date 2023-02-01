///////////////////////////////////////////////////////////////////////////////////////
/// \file getcliminput.cpp
/// \brief Input module for driver file generated by GetClim tool
///
/// \author Ben Smith
/// $Date: 2016-07-26 16:25:45 +0100 (Tue, 26 Jul 2016) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "getcliminput.h"
#include "soilinput.h"
#include "guessmath.h"

#include "driver.h"
#include "outputchannel.h"
#include <stdio.h>

REGISTER_INPUT_MODULE("getclim", GetclimInput)

// Anonymous namespace for variables and functions with file scope
namespace {

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double* mtemp, double* mprec, double* msun, double* mdtr,
					double* dtemp, double* dprec, double* dsun, double* ddtr) {
	interp_monthly_means_conserve(mtemp, dtemp);
	interp_monthly_totals_conserve(mprec, dprec, 0);
	interp_monthly_means_conserve(msun, dsun, 0, 370);
	interp_monthly_means_conserve(mdtr, ddtr, 0);
}

} // namespace

GetclimInput::GetclimInput() {

	in_clim = NULL;
}


// Reads header of GetClim driver file
void GetclimInput::init_climate(double& dlon, double& dlat)
{
	xtring text = "";
	int i = 0;
	//Dummy for reading non-used variables from climate file
	int nyear_transient,nyear_stabilisation;
	double scenario_delta_temp, scenario_factor_precip, scenario_factor_rad, scenario_factor_co2, scenario_factor_Ndep;

	readfor(in_clim, ""); // Header
	readfor(in_clim, ""); // Blank line
	readfor(in_clim, ""); // Header2
	readfor(in_clim, ""); // Blank line
	readfor(in_clim, ""); // Column headings
	readfor(in_clim, "f,f,i,i,i,i,f,f,f,f", &dlon, &dlat, &soilcode, &nyear_spinup, &nyear_transient, &nyear_stabilisation,
		&scenario_delta_temp, &scenario_factor_precip, &scenario_factor_rad, &scenario_factor_co2, &scenario_factor_Ndep);
	readfor(in_clim, ""); // Blank line
	readfor(in_clim, ""); // Column headings
	readfor(in_clim, ""); // Column headings2

	nyear_scenario = nyear_transient + nyear_stabilisation;

	if (ifbvoc)
		dprintf("WARNING: No data available for dtr in sample data set!\nNo daytime temperature correction for BVOC calculations applied.");

	if (nyear_spinup <= freenyears) {
		fail("Error: freenyears must be smaller than nyear_spinup");
	}
}	/////////////////////////////////////////////////////////////////////

void GetclimInput::init() {

	// DESCRIPTION
	// Initialises input (e.g. opening files)

	// Getclim input module currently only works with the old INTERP weather generator and GLOBFIRM (or NO FIRE).
	if (weathergenerator == GWGEN || firemodel == BLAZE) {
		fail("Getclim input module currently only works with the INTERP weather generator and the fire model GLOBFIRM (or no fire with NOFIRE).\n Make sure that both of them are set correctly in global.ins.");
	}

	// Open landcover files. May reduce pftlist, stlist and mtlist. Must be called before management_input->init()
	landcover_input.init();
	// Open management files
	management_input.init();
	// Open additional files
	misc_input.init();

	xtring driver_file_path = param["getclim_driver_file"].str;

	// Open GetClim driver file
	in_clim = fopen(driver_file_path, "rt");
	if (!in_clim) fail("Could not open %s for input\n", (char*)driver_file_path);
	init_climate(lon, lat);

	// Initialise grid cell count (only one grid cell per simulation)
	grid_count = 0;

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

bool GetclimInput::getgridcell(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	bool LUerror=false;

	// Load land use data for this grid cell from files

	if (grid_count++ > 0) return false;
			
	if(readdisturbance || readdisturbance_st || readelevation_st) {
		// Not all gridcells have to be included in input file
		misc_input.loaddisturbance(lon, lat);
		misc_input.loadelevation(lon, lat);
	}

	if(run_landcover) {
		LUerror = landcover_input.loadlandcover(lon, lat);
		if(!LUerror)
			LUerror = management_input.loadmanagement(lon, lat);
	}
	
	if (LUerror) return false;

	// else ...

	dprintf("\nCommencing simulation for stand at (%g,%g)\n\n",lon,lat);
	// Tell framework the coordinates of this grid cell
	gridcell.set_coordinates(lon, lat);
	// The insolation data will be sent (in function getclimate, below)
	// as total shortwave radiation for whole time step (based on CRU-NCEP database)
	gridcell.climate.instype=SWRAD_TS;

	// Tell framework the soil type of this grid cell
	soil_parameters(gridcell.soiltype,soilcode);

	// For Windows shell - clear graphical output
	clear_all_graphs();

	return true; // simulate this stand
}


void GetclimInput::getlandcover(Gridcell& gridcell) {

	landcover_input.getlandcover(gridcell);
	landcover_input.get_land_transitions(gridcell);
}

bool GetclimInput::getclimate(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	double progress;
	double mprec[12], mtemp[12], msun[12],mwet[12]; 
	double mdtr[12] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	int dummy;
	xtring text;

	Climate& climate = gridcell.climate;

	// Retrieve data for one year on first day of year

	if (date.day == 0) {

		// Return false if last year was the last for the simulation
		if (date.year == nyear_spinup + nyear_scenario) {

			dprintf("\nModel run terminated\n");
			return false;
		}

		readfor(in_clim, "a,i,f,f,12f,12f,12f,12f", &text, &dummy, &co2, &ndep, mtemp, mprec, msun, mwet);

		// Interpolate monthly values to daily values
		interp_climate(mtemp, mprec, msun, mdtr, dtemp, dprec, dsun, ddtr);

		// Recalculate precipitation values using weather generator
		if (ifrainonwetdaysonly) {
			prdaily(mprec, dprec, mwet, gridcell.seed);
		}

	}

	// Send environmental values for today to framework

	gridcell.dNH4dep  = ndep / 2.0 / 365.0 * HA_PER_M2;
	gridcell.dNO3dep = ndep / 2.0 / 365.0 * HA_PER_M2;

	climate.co2 = co2;

	climate.temp  = dtemp[date.day];
	climate.prec  = dprec[date.day];
	climate.insol = dsun[date.day];

	// Daily temperature range (placeholder code, required by BVOC calculations, not available in GetClim driver file)

	climate.dtr=ddtr[date.day];
	
	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(date.year)/(double)(nyear_spinup+nyear_scenario);

			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}

GetclimInput::~GetclimInput() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

	fclose(in_clim);
}
