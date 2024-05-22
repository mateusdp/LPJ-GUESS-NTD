////////////////////////////////////////////////////////////////////////////////////////
/// \file cfinput.cpp
/// \brief Input module for CF conforming NetCDF files
///
/// \author Joe Siltberg
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

// Include N-deposition data read in as netCDF (author: Tim Anders, with adaptations by Mateus Dantas)
#include "config.h"
#include "cfinput.h"

#ifdef HAVE_NETCDF

#include "guess.h"
#include "driver.h"
#include "weathergen.h"
#include "guessstring.h"
#include <fstream>
#include <sstream>
#include <algorithm>

REGISTER_INPUT_MODULE("cf", CFInput)

using namespace GuessNC::CF;

namespace {

const int SECONDS_PER_DAY = 24*60*60;

// Converts a CF standard name to one of our insolation types
// Calls fail() if the standard name is invalid
insoltype cf_standard_name_to_insoltype(const std::string& standard_name) {
	if (standard_name == "surface_downwelling_shortwave_flux_in_air" ||
	    standard_name == "surface_downwelling_shortwave_flux") {
		return SWRAD_TS;
	}
	else if (standard_name == "surface_net_downward_shortwave_flux") {
		return NETSWRAD_TS;
	}
	else if (standard_name == "cloud_area_fraction") {
		return SUNSHINE;
	}
	else {
		fail("Unknown insolation type: %s", standard_name.c_str());
		return SUNSHINE; // To avoid compiler warning
	}
}

// Gives the maximum allowed value for insolation, given an insolation type
// Used as an upper limit when interpolating from monthly to daily values
double max_insolation(insoltype instype) {
	if (instype == SUNSHINE) {
		return 100;
	}
	else {
		return std::numeric_limits<double>::max();
	}
}

// Checks if a DateTime is at the first day of the year
bool first_day_of_year(GuessNC::CF::DateTime dt) {
	return dt.get_month() == 1 && dt.get_day() == 1;
}

// Checks if a DateTime is in January
bool first_month_of_year(GuessNC::CF::DateTime dt) {
	return dt.get_month() == 1;
}

// Compares a Date with a GuessNC::CF::DateTime to see if the Date is on an earlier day
bool earlier_day(const Date& date, int calendar_year,
                 const GuessNC::CF::DateTime& date_time) {
	std::vector<int> d1(3),d2(3);

	d1[0] = calendar_year;
	d2[0] = date_time.get_year();

	d1[1] = date.month+1;
	d2[1] = date_time.get_month();

	d1[2] = date.dayofmonth+1;
	d2[2] = date_time.get_day();

	return d1 < d2;
}

// Compares a Date with a GuessNC::CF::DateTime to see if the Date is on a later day
// The date object must know about its calendar years (i.e. set_first_calendar_year must
// have been called)
bool later_day(const Date& date,
               const GuessNC::CF::DateTime& date_time) {
	std::vector<int> d1(3),d2(3);

	d1[0] = date.get_calendar_year();
	d2[0] = date_time.get_year();

	d1[1] = date.month+1;
	d2[1] = date_time.get_month();

	d1[2] = date.dayofmonth+1;
	d2[2] = date_time.get_day();

	return d1 > d2;
}

// Checks if the variable contains daily data
bool is_daily(const GuessNC::CF::GridcellOrderedVariable* cf_var) {

	// Check if first and second timestep is one day apart

	DateTime dt1 = cf_var->get_date_time(0);
	DateTime dt2 = cf_var->get_date_time(1);

	dt1.add_time(1, GuessNC::CF::DAYS, cf_var->get_calendar_type());

	return dt1 == dt2;
}

// Returns a DateTime in the last day for which the variable has data.
// For daily data, this is simply the day of the last timestep, for monthly data
// we need to find the last day of the last timestep's month.
GuessNC::CF::DateTime last_day_to_simulate(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	GuessNC::CF::DateTime last = cf_var->get_date_time(cf_var->get_timesteps()-1);
	if (is_daily(cf_var)) {
		return last;
	}
	else {
		// Not daily, assume monthly.
		GuessNC::CF::DateTime prev = last;
		GuessNC::CF::DateTime next = last;

		do {
			prev = next;
			next.add_time(1, GuessNC::CF::DAYS, cf_var->get_calendar_type());
		} while (next.get_month() == last.get_month());

		return prev;
	}
}

// Verifies that a CF variable with air temperature data contains what we expect
void check_temp_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	if (cf_var->get_standard_name() != "air_temperature") {
		fail("Temperature variable doesn't seem to contain air temperature data");
	}
	if (cf_var->get_units() != "K") {
		fail("Temperature variable doesn't seem to be in Kelvin");
	}
}

// Verifies that a CF variable with precipitation data contains what we expect
void check_prec_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	if (cf_var->get_standard_name() == "precipitation_flux") {
		if (cf_var->get_units() != "kg m-2 s-1") {
			fail("Precipitation is given as flux but does not have the correct unit (kg m-2 s-1)");
		}
	}
	else if (cf_var->get_standard_name() == "precipitation_amount") {
		if (cf_var->get_units() != "kg m-2") {
			fail("Precipitation is given as amount but does not have the correct unit (kg m-2)");
		}
	}
	else {
		fail("Unrecognized precipitation type");
	}
}

// Verifies that a CF variable with insolation data contains what we expect
void check_insol_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	if (cf_var->get_standard_name() != "surface_downwelling_shortwave_flux_in_air" &&
	    cf_var->get_standard_name() != "surface_downwelling_shortwave_flux" &&
	    cf_var->get_standard_name() != "surface_net_downward_shortwave_flux" &&
	    cf_var->get_standard_name() != "cloud_area_fraction") {
		fail("Insolation variable doesn't seem to contain insolation data");
	}

	if (cf_var->get_standard_name() == "cloud_area_fraction") {
		if (cf_var->get_units() != "1") {
			fail("Unrecognized unit for cloud cover");
		}
	}
	else {
		if (cf_var->get_units() != "W m-2") {
			fail("Insolation variable given as radiation but unit doesn't seem to be in W m-2");
		}
	}
}

// Verifies that a CF variable with wetdays data contains what we expect
void check_wetdays_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* wetdays_standard_name =
		"number_of_days_with_lwe_thickness_of_precipitation_amount_above_threshold";

	if (cf_var && cf_var->get_standard_name() != wetdays_standard_name) {
		fail("Wetdays variable should have standard name %s", wetdays_standard_name);
	}
}

// Verifies that a CF variable with pressure data contains what we expect
void check_pres_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* pres_standard_name = "surface_air_pressure";
	if (cf_var->get_standard_name() != pres_standard_name) {
		fail("Pressure variable should have standard name %s ",pres_standard_name);
	}
	if (cf_var->get_units() != "Pa") {
		fail("Pressure must be given in Pa!");
	}
}

// Verifies that a CF variable with specific humidity data contains what we expect
void check_specifichum_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* standard_name = "specific_humidity";
	if (cf_var->get_standard_name() != standard_name) {
		fail("QAir variable should have standard name %s ",standard_name);
	}
	if (cf_var->get_units() != "1") {
		fail("Specific Humidity must be dimensionless (here, '1'!");
	}
}

// Verifies that a CF variable with relative humidity data contains what we expect
void check_relhum_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* standard_name = "relative_humidity";
	if (cf_var->get_standard_name() != standard_name) {
		fail("Relative humidity variable should have standard name %s ",standard_name);
	}
	if (cf_var->get_units() != "1") {
		fail("Relative Humidity must be dimensionless (here, '1'!");
	}
}

// Verifies that a CF variable with wind-speed data contains what we expect
void check_wind_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* standard_name = "wind_speed";
	if (cf_var->get_standard_name() != standard_name) {
		fail("Wind variable should have standard name %s ",standard_name);
	}
	if (cf_var->get_units() != "m s-1") {
		fail("Wind data must be given in m s-1 !");
	}
}

// Verifies that a CF variable with dry NHx deposition data contains what we expect
void check_NHxdrydep_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* standard_name = "mNHxdrydep";
	if (cf_var->get_standard_name() != standard_name) {
		fail("Dry NHx deposition variable should have standard name %s ", standard_name);
	}
	if (cf_var->get_units() != "kg N m-2 d-1") { //Changed unit by Dantas, from s-1 to d-1
		fail("Dry NHx deposition data must be given in kg N m-2 d-1 !");  //Changed unit by Dantas
	}
	else {
		dprintf("All fine with NHxdrydep netCDF file.\n");
	}
} //TA added check for NHxdrydep data (Changed unit by Dantas)


// Checks if two variables contain data for the same time period
//
// Compares start and end of time series, the day numbers are only compared if
// both variables are daily.
void check_compatible_timeseries(const GuessNC::CF::GridcellOrderedVariable* var1,
                                 const GuessNC::CF::GridcellOrderedVariable* var2) {
	GuessNC::CF::DateTime start1, start2, end1, end2;

	const std::string error_message = format_string("%s and %s have incompatible timeseries",
		var1->get_variable_name().c_str(), var2->get_variable_name().c_str());

	start1 = var1->get_date_time(0);
	start2 = var2->get_date_time(0);

	end1 = var1->get_date_time(var1->get_timesteps() - 1);
	end2 = var2->get_date_time(var2->get_timesteps() - 1);

	if (start1.get_year() != start2.get_year() ||
		start1.get_month() != start2.get_month()) {
		fail(error_message.c_str());
	}

	if (end1.get_year() != end2.get_year() ||
		end1.get_month() != end2.get_month()) {
		fail(error_message.c_str());
	}

	if (is_daily(var1) && is_daily(var2)) {
		if (start1.get_day() != start2.get_day() ||
			end1.get_day() != end2.get_day()) {
			fail(error_message.c_str());
		}
	}
}

// Makes sure all variables have compatible time series
void check_compatible_timeseries(const std::vector<GuessNC::CF::GridcellOrderedVariable*> variables) {

	for (size_t i = 0; i < variables.size(); ++i) {
		for (size_t j = i + 1; j < variables.size(); ++j) {
			check_compatible_timeseries(variables[i], variables[j]);
		}
	}
}

void check_same_spatial_domains(const std::vector<GuessNC::CF::GridcellOrderedVariable*> variables) {

	for (size_t i = 1; i < variables.size(); ++i) {
		if (!variables[0]->same_spatial_domain(*variables[i])) {
			fail("%s and %s don't have the same spatial domain",
				variables[0]->get_variable_name().c_str(),
				variables[i]->get_variable_name().c_str());
		}
	}
}

// Compute relative humidity from specific humidity, temperature and pressure
double calc_relative_humidity(double temp, double specific_humidity, double pressure) {

	// qair  specific humidity, dimensionless (e.g. kg/kg) 
	// temp  temperature in degrees C
	// press pressure in Pa
	// rh    relative humidity in frac.
	if ( pressure > 106000 || pressure < 10000 ) {
		fail("Unit for pressure must be [Pa]: calc_relative_humidity(cfinput.cpp)");
	} 
	if ( temp  > 80. ) {
		fail("Unit for temperature must be [deg C]: calc_relative_humidity(cfinput.cpp)");
	} 
	double pres_hPa = pressure / 100.; // convert to hPa

	// saturation water-vapour pressure following August-Roche-Magnus Formula
	double es   = 6.112 * exp(17.67 * temp/(temp + 243.5));

	// water-vapour pressure
	// derived from approximation for s = rho_w/(rho_dryAir - rho_w)    
	double e    = specific_humidity * pres_hPa / (0.378 * specific_humidity + 0.622);
	double rh   = min(max(e / es ,0.),1.) ;
	return rh;		
}

}

CFInput::CFInput()
	: cf_temp(0),
	  cf_prec(0),
	  cf_insol(0),
	  cf_wetdays(0),
	  cf_min_temp(0),
	  cf_max_temp(0),
	  cf_pres(0),
	  cf_specifichum(0),
	  cf_relhum(0),
	  cf_wind(0),
	  cf_mNHxdrydep(0),
	  cf_mNOydrydep(0),
	  cf_mNHxwetdep(0),
	  cf_mNOywetdep(0),
	  cf_mPOxdrydep(0),
	  cf_mPOxwetdep(0),
	  ndep_timeseries("historic") {

	// Declare instruction file parameters
	declare_parameter("ndep_timeseries", &ndep_timeseries, 10, "Nitrogen deposition time series to use (historic, rcp26, rcp45, rcp60 or rcp85");
	 SoilInput soilinput;
}

CFInput::~CFInput() {
	delete cf_temp;
	delete cf_prec;
	delete cf_insol;
	delete cf_wetdays;
	delete cf_min_temp;
	delete cf_max_temp;
	delete cf_pres;
	delete cf_specifichum;
	delete cf_relhum;
	delete cf_wind;
	delete cf_mNHxdrydep;
	delete cf_mNOydrydep;
	delete cf_mNHxwetdep;
	delete cf_mNOywetdep;
	delete cf_mPOxdrydep;
	delete cf_mPOxwetdep;


	cf_temp = 0;
	cf_prec = 0;
	cf_insol = 0;
	cf_wetdays = 0;
	cf_min_temp = 0;
	cf_max_temp = 0;
	cf_pres = 0;
	cf_specifichum = 0;
	cf_relhum = 0;
	cf_wind = 0;
	cf_mNHxdrydep = 0;
	cf_mNOydrydep = 0;
	cf_mNHxwetdep = 0;
	cf_mNOywetdep = 0;
	cf_mPOxdrydep = 0;
	cf_mPOxwetdep = 0;
}

void CFInput::init() {

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	// Try to open the NetCDF files
	try {
		cf_temp = new GridcellOrderedVariable(param["file_temp"].str, param["variable_temp"].str);
		cf_prec = new GridcellOrderedVariable(param["file_prec"].str, param["variable_prec"].str);
		cf_insol = new GridcellOrderedVariable(param["file_insol"].str, param["variable_insol"].str);

		if (param["file_wetdays"].str != "") {
			cf_wetdays = new GridcellOrderedVariable(param["file_wetdays"].str, param["variable_wetdays"].str);
		}

		if (param["file_min_temp"].str != "") {
			cf_min_temp = new GridcellOrderedVariable(param["file_min_temp"].str, param["variable_min_temp"].str);
		}

		if (param["file_max_temp"].str != "") {
			cf_max_temp = new GridcellOrderedVariable(param["file_max_temp"].str, param["variable_max_temp"].str);
		}

		if (param["file_pres"].str != "") {
			cf_pres = new GridcellOrderedVariable(param["file_pres"].str, param["variable_pres"].str);
		}

		if (param["file_specifichum"].str != "") {
			cf_specifichum = new GridcellOrderedVariable(param["file_specifichum"].str, param["variable_specifichum"].str);
		}

		if (param["file_relhum"].str != "") {
			cf_relhum = new GridcellOrderedVariable(param["file_relhum"].str, param["variable_relhum"].str);
		}

		if (param["file_wind"].str != "") {
			cf_wind = new GridcellOrderedVariable(param["file_wind"].str, param["variable_wind"].str);
		}

		if (param.isparam("file_mNHxdrydep") && param["file_mNHxdrydep"].str != "") {
			cf_mNHxdrydep = new GridcellOrderedVariable(param["file_mNHxdrydep"].str, param["variable_mNHxdrydep"].str);
			dprintf("Opening of mNHxdrydep netcdf successful.\n"); //TA test opening of NHxdrydep netcdf data
		}

		if (param.isparam("file_mNOydrydep") && param["file_mNOydrydep"].str != "") {
			cf_mNOydrydep = new GridcellOrderedVariable(param["file_mNOydrydep"].str, param["variable_mNOydrydep"].str);
		}

		if (param.isparam("file_mNHxwetdep") && param["file_mNHxwetdep"].str != "") {
			cf_mNHxwetdep = new GridcellOrderedVariable(param["file_mNHxwetdep"].str, param["variable_mNHxwetdep"].str);
		}

		if (param.isparam("file_mNOywetdep") && param["file_mNOywetdep"].str != "") {
			cf_mNOywetdep = new GridcellOrderedVariable(param["file_mNOywetdep"].str, param["variable_mNOywetdep"].str);
		}

		if (param.isparam("file_mPOxdrydep") && param["file_mPOxdrydep"].str != "") {
			cf_mPOxdrydep = new GridcellOrderedVariable(param["file_mPOxdrydep"].str, param["variable_mPOxdrydep"].str);
		}

		if (param.isparam("file_mPOxwetdep") && param["file_mPOxwetdep"].str != "") {
			cf_mPOxwetdep = new GridcellOrderedVariable(param["file_mPOxwetdep"].str, param["variable_mPOxwetdep"].str);
		}

	}
	catch (const std::runtime_error& e) {
		fail(e.what());
	}

	// Make sure they contain what we expect

	check_temp_variable(cf_temp);

	check_prec_variable(cf_prec);

	check_insol_variable(cf_insol);

	check_wetdays_variable(cf_wetdays);

	if (cf_min_temp) {
		check_temp_variable(cf_min_temp);
	}

	if (cf_max_temp) {
		check_temp_variable(cf_max_temp);
	}

	if (cf_pres) {
		check_pres_variable(cf_pres);
	}

	if (cf_specifichum) {
		check_specifichum_variable(cf_specifichum);
	}

	if (cf_relhum) {
		check_relhum_variable(cf_relhum);
	}

	if (cf_wind) {
		check_wind_variable(cf_wind);
	}

	check_compatible_timeseries(all_variables());

	check_same_spatial_domains(all_variables());

	if (cf_mNHxdrydep) {
		if (!cf_mNHxdrydep->same_spatial_domain(*cf_temp)) {
			fail("dry ndep %s has not the same spatial domain as climate data", (char*)param["file_mNHxdrydep"].str);
		};
	}

	if (cf_mNOydrydep) {
		if (!cf_mNOydrydep->same_spatial_domain(*cf_temp)) {
			fail("dry ndep %s has not the same spatial domain as climate data", (char*)param["file_mNOydrydep"].str);
		};
	}

	if (cf_mNHxwetdep) {
		if (!cf_mNHxwetdep->same_spatial_domain(*cf_temp)) {
			fail("wet ndep %s has not the same spatial domain as climate data", (char*)param["file_mNHxwetdep"].str);
		};
	}

	if (cf_mNOywetdep) {
		if (!cf_mNOywetdep->same_spatial_domain(*cf_temp)) {
			fail("wet ndep %s has not the same spatial domain as climate data", (char*)param["file_mNOywetdep"].str);
		};
	}

	if (cf_mPOxdrydep) {
		if (!cf_mPOxdrydep->same_spatial_domain(*cf_temp)) {
			fail("dry pdep %s has not the same spatial domain as climate data", (char*)param["file_mPOxdrydep"].str);
		};
	}
	if (cf_mPOxwetdep) {
		if (!cf_mPOxwetdep->same_spatial_domain(*cf_temp)) {
			fail("wet pdep %s has not the same spatial domain as climate data", (char*)param["file_mPOxwetdep"].str);
		};
	}

	extensive_precipitation = cf_prec->get_standard_name() == "precipitation_amount";

	// Read list of localities and store in gridlist member variable

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist_cf"].str;

	std::ifstream ifs(file_gridlist, std::ifstream::in);

	if (!ifs.good()) fail("CFInput::init: could not open %s for input",(char*)file_gridlist);

	std::string line;
	while (getline(ifs, line)) {

		// Read next record in file
		int rlat, rlon;
		int landid;
		std::string descrip;
		Coord c;

		std::istringstream iss(line);

		if (cf_temp->is_reduced()) {
			if (iss >> landid) {
				getline(iss, descrip);

				c.landid = landid;
			}
		}
		else {
			if (iss >> rlon >> rlat) {
				getline(iss, descrip);

				c.rlat = rlat;
				c.rlon = rlon;

			}
			else {
				fail("The gridlist for netCDF input must be in X,Y coordinates");
			}
		}
		c.descrip = (xtring)trim(descrip).c_str();
		gridlist.push_back(c);
	}

	current_gridcell = gridlist.begin();

	// Open landcover files. May reduce pftlist, stlist and mtlist. Must be called before management_input->init()
	landcover_input.init();
	// Open management files
	management_input.init();
	// Open additional files
	misc_input.init();

	date.set_first_calendar_year(cf_temp->get_date_time(0).get_year() - nyear_spinup);


	soilinput.init(param["file_soildata"].str);

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

bool CFInput::getgridcell(Gridcell& gridcell) {

	double lon, lat;
	double cru_lon, cru_lat;
	int soilstatus = 1;

	// Load data for next gridcell, or if that fails, skip ahead until
	// we find one that works.
	while (current_gridcell != gridlist.end() &&
	       !load_data_from_files(lon, lat)){
		++current_gridcell;
	}

	/*cru_lon = floor(lon * 2.0) / 2.0 + 0.25;
	cru_lat = floor(lat * 2.0) / 2.0 + 0.25;*/

	if (current_gridcell == gridlist.end()) {
		// simulation finished
		return false;
	}

	cru_lon = floor(lon * 2.0) / 2.0 + 0.25;
	cru_lat = floor(lat * 2.0) / 2.0 + 0.25;

	if(readdisturbance || readdisturbance_st || readelevation_st) {
		// Not all gridcells have to be included in input file
		misc_input.loaddisturbance(lon, lat);
		misc_input.loadelevation(lon, lat);
	}

//	gridcell.climate.mean_elevation = elevation;		// Get elevation from cru_ncep
//	if(readelevation_st)
//		dprintf("Mean elevation = %d\n", elevation);

	if (run_landcover) {
		bool LUerror = false;
		LUerror = landcover_input.loadlandcover(lon, lat);
		if (!LUerror)
			LUerror = management_input.loadmanagement(lon, lat);
		if (LUerror) {
			dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n", lon, lat);
			return false;
		}
	}

	gridcell.set_coordinates(lon, lat);

	// Load spinup data for all variables

	load_spinup_data(cf_temp, spinup_temp);
	load_spinup_data(cf_prec, spinup_prec);
	load_spinup_data(cf_insol, spinup_insol);

	if (cf_wetdays) {
		load_spinup_data(cf_wetdays, spinup_wetdays);
	}

	if (cf_min_temp) {
		load_spinup_data(cf_min_temp, spinup_min_temp);
	}

	if (cf_max_temp) {
		load_spinup_data(cf_max_temp, spinup_max_temp);
	}

	if (cf_pres) {
		load_spinup_data(cf_pres, spinup_pres);
	}

	if (cf_specifichum) {
		load_spinup_data(cf_specifichum, spinup_specifichum);
	}

	if (cf_relhum) {
		load_spinup_data(cf_relhum, spinup_relhum);
	}

	if (cf_wind) {
		load_spinup_data(cf_wind, spinup_wind);
	}
	
	spinup_temp.detrend_data();

	gridcell.climate.instype = cf_standard_name_to_insoltype(cf_insol->get_standard_name());

	// Get nitrogen deposition, using the found CRU coordinates
	/* Since the historic data set does not reach decade 2010-2019,
	* we need to use the RCP data for the last decade. */
	//ndep.getndep(param["file_ndep"].str, cru_lon, cru_lat, Lamarque::RCP60);

	if (param["file_ndep"].str != "") {
		double cru_lon;
		double cru_lat;
		// this tries to calculate the 0.5deg cru lon lat, might not work if that lon lat is not within the ndep data
		cru_lon = floor(lon * 2.0) / 2.0 + 0.25;
		cru_lat = floor(lat * 2.0) / 2.0 + 0.25;
		//dprintf("lon = %3.2f, lat = %3.2f)\n", lon, lat); // TA inserted print
		//dprintf("Cru_lon = %3.2f, Cru_lat = %3.2f)\n", cru_lon, cru_lat); // TA inserted print

		dprintf("Using Nitrogen deposition for (%3.2f,%3.2f)\n", cru_lon, cru_lat);
		// Get nitrogen deposition, using the estimated CRU coordinates
		//ndep.getndep(param["file_ndep"].str, cru_lon, cru_lat,Lamarque::parse_timeseries(ndep_timeseries));
		ndep.getndep(param["file_ndep"].str, cru_lon, cru_lat, Lamarque::RCP60); //TA parse_timeseries is not working, have to use RCP60
	}


	soilinput.get_soil(lon, lat, gridcell);

	historic_timestep_temp = -1;
	historic_timestep_prec = -1;
	historic_timestep_insol = -1;
	historic_timestep_wetdays = -1;
	historic_timestep_min_temp = -1;
	historic_timestep_max_temp = -1;

	historic_timestep_pres = -1;
	historic_timestep_specifichum = -1;
	historic_timestep_relhum = -1;
	historic_timestep_wind = -1;

	dprintf("\nCommencing simulation for gridcell at (%g,%g)\n", lon, lat);
	if (current_gridcell->descrip != "") {
		dprintf("Description: %s\n", current_gridcell->descrip.c_str());
	}
	dprintf("Using soil code and Nitrogen deposition for (%3.1f,%3.1f)\n", cru_lon, cru_lat);

	return true;
}

bool CFInput::load_data_from_files(double& lon, double& lat){

	int rlon = current_gridcell->rlon;
	int rlat = current_gridcell->rlat;
	int landid = current_gridcell->landid;

	// Try to load the data from the NetCDF files

	if (cf_temp->is_reduced()) {
		if (!cf_temp->load_data_for(landid) ||
		    !cf_prec->load_data_for(landid) ||
		    !cf_insol->load_data_for(landid) ||
		    (cf_wetdays && !cf_wetdays->load_data_for(landid)) ||
		    (cf_min_temp && !cf_min_temp->load_data_for(landid)) ||
		    (cf_max_temp && !cf_max_temp->load_data_for(landid)) ||
		    (cf_pres && !cf_pres->load_data_for(landid)) ||
		    (cf_specifichum && !cf_specifichum->load_data_for(landid)) ||
		    (cf_relhum && !cf_relhum->load_data_for(landid)) ||
		    (cf_wind && !cf_wind->load_data_for(landid)) ||
			(cf_mNHxdrydep && true) || //TA landid not implemented for the ff. nc's, hence fail
			(cf_mNOydrydep && true) ||
			(cf_mNHxwetdep && true) ||
			(cf_mNOywetdep && true) ||
			(cf_mPOxdrydep && true) ||
			(cf_mPOxwetdep && true)){
			dprintf("Failed to load data for (%d) from NetCDF files, skipping.\n", landid);
			return false;
		}
	}
	else {
		if (!cf_temp->load_data_for(rlon, rlat) ||
		    !cf_prec->load_data_for(rlon, rlat) ||
		    !cf_insol->load_data_for(rlon, rlat) ||
		    (cf_wetdays && !cf_wetdays->load_data_for(rlon, rlat)) ||
		    (cf_min_temp && !cf_min_temp->load_data_for(rlon, rlat)) ||
		    (cf_max_temp && !cf_max_temp->load_data_for(rlon, rlat))||
		    (cf_pres && !cf_pres->load_data_for(rlon, rlat))||
		    (cf_specifichum && !cf_specifichum->load_data_for(rlon, rlat))||
		    (cf_relhum && !cf_relhum->load_data_for(rlon, rlat))||
		    (cf_wind && !cf_wind->load_data_for(rlon, rlat)) ||
			(cf_mNHxdrydep && !cf_mNHxdrydep->load_data_for(rlon, rlat)) ||
			(cf_mNOydrydep && !cf_mNOydrydep->load_data_for(rlon, rlat)) ||
			(cf_mNHxwetdep && !cf_mNHxwetdep->load_data_for(rlon, rlat)) ||
			(cf_mNOywetdep && !cf_mNOywetdep->load_data_for(rlon, rlat)) ||
			(cf_mPOxdrydep && !cf_mPOxdrydep->load_data_for(rlon, rlat)) || 
			(cf_mPOxwetdep && !cf_mPOxwetdep->load_data_for(rlon, rlat))) {
			dprintf("Failed to load data for (%d, %d) from NetCDF files, skipping.\n", rlon, rlat);
			return false;
		}
	}

	// Get lon/lat for the gridcell

	if (cf_temp->is_reduced()) {
		cf_temp->get_coords_for(landid, lon, lat);
	}
	else {
		cf_temp->get_coords_for(rlon, rlat, lon, lat);
	}

	return true;
}

void CFInput::get_yearly_data(std::vector<double>& data,
                              const GenericSpinupData& spinup,
                              GridcellOrderedVariable* cf_historic,
                              int& historic_timestep) {
	// Extract all values for this year, for one variable,
	// either from spinup dataset or historical dataset

	int calendar_year = date.get_calendar_year();

	if (is_daily(cf_historic)) {

		data.resize(date.year_length());

		// This function is called at the first day of the year, so current_day
		// starts at Jan 1, then we step through the whole year, getting data
		// either from spinup or historical period.
		Date current_day = date;

		while (current_day.year == date.year) {

			// In the spinup?
			if (earlier_day(current_day, calendar_year, cf_historic->get_date_time(0))) {

				int spinup_day = current_day.day;
				// spinup object always has 365 days, deal with leap years
				if (current_day.ndaymonth[1] == 29 && current_day.month > 1) {
					--spinup_day;
				}
				data[current_day.day]  = spinup[spinup_day];
			}
			else {
				// Historical period

				if (historic_timestep + 1 < cf_historic->get_timesteps()) {

					++historic_timestep;
					GuessNC::CF::DateTime dt = cf_historic->get_date_time(historic_timestep);

					// Deal with calendar mismatch

					// Leap day in NetCDF variable but not in LPJ-GUESS?
					if (dt.get_month() == 2 && dt.get_day() == 29 &&
						current_day.ndaymonth[1] == 28) {
						++historic_timestep;
					}
					// Leap day in LPJ-GUESS but not in NetCDF variable?
					else if (current_day.month == 1 && current_day.dayofmonth == 28 &&
						cf_historic->get_calendar_type() == NO_LEAP) {
						--historic_timestep;
					}
				}

				if (historic_timestep < cf_historic->get_timesteps()) {
					data[current_day.day]  = cf_historic->get_value(max(0, historic_timestep));
				}
				else {
					// Past the end of the historical period, these days wont be simulated.
					data[current_day.day] = data[max(0, current_day.day-1)];
				}
			}

			current_day.next();
		}
	}
	else {

		// for now, assume that data set must be monthly since it isn't daily

		data.resize(12);

		for (int m = 0; m < 12; ++m) {

			GuessNC::CF::DateTime first_date = cf_historic->get_date_time(0);

			// In the spinup?
			if (calendar_year < first_date.get_year() ||
			    (calendar_year == first_date.get_year() &&
			     m+1 < first_date.get_month())) {
				data[m] = spinup[m];
			}
			else {
				// Historical period
				if (historic_timestep + 1 < cf_historic->get_timesteps()) {
					++historic_timestep;
				}

				if (historic_timestep < cf_historic->get_timesteps()) {
					data[m] = cf_historic->get_value(historic_timestep);
				}
				else {
					// Past the end of the historical period, these months wont be simulated.
					data[m] = data[max(0, m-1)];
				}
			}
		}
	}
}

void CFInput::populate_daily_array(double* daily,
                                   const GenericSpinupData& spinup,
                                   GridcellOrderedVariable* cf_historic,
                                   int& historic_timestep,
                                   double minimum,
                                   double maximum) {

	// Get the data from spinup and/or historic
	std::vector<double> data;
	get_yearly_data(data, spinup, cf_historic, historic_timestep);

	if (is_daily(cf_historic)) {
		// Simply copy from data to daily

		std::copy(data.begin(), data.end(), daily);
	}
	else {
		// for now, assume that data set must be monthly since it isn't daily

		// Interpolate from monthly to daily values

		interp_monthly_means_conserve(&data.front(), daily, minimum, maximum);
	}
}

void CFInput::populate_daily_prec_array(long& seed) {

	// Get the data from spinup and/or historic
	std::vector<double> prec_data;
	get_yearly_data(prec_data, spinup_prec, cf_prec, historic_timestep_prec);

	std::vector<double> wetdays_data;
	if (cf_wetdays) {
		get_yearly_data(wetdays_data, spinup_wetdays, cf_wetdays, historic_timestep_wetdays);
	}

	if (is_daily(cf_prec)) {
		// Simply copy from data to daily, and if needed convert from
		// precipitation rate to precipitation amount

		for (size_t i = 0; i < prec_data.size(); ++i) {
			dprec[i] = prec_data[i];

			if (!extensive_precipitation) {
				dprec[i] *= SECONDS_PER_DAY;
			}
		}
	}
	else {
		// for now, assume that data set must be monthly since it isn't daily

		// If needed convert from precipitation rate to precipitation amount
		if (!extensive_precipitation) {
			for (int m = 0; m < 12; ++m) {
				// TODO: use the dataset's calendar type to figure out number of days in month?
				prec_data[m] *= SECONDS_PER_DAY * date.ndaymonth[m];
			}
		}

		if (cf_wetdays) {
			prdaily(&prec_data.front(), dprec, &wetdays_data.front(), seed);
		}
		else {
			interp_monthly_totals_conserve(&prec_data.front(), dprec, 0);
		}
	}
}

void CFInput::populate_daily_arrays(Gridcell& gridcell) {
	// Extract daily values for all days in this year, either from
	// spinup dataset or historical dataset

	if ( !is_daily(cf_temp) && weathergenerator == GWGEN ) {

		int instype = cf_standard_name_to_insoltype(cf_insol->get_standard_name());
		
		// TODO IMPLEMENT cloud-frac
		if (!cf_min_temp || !cf_max_temp || !cf_wind || ( ( !cf_pres || !cf_specifichum ) && !cf_relhum ) ||
		    instype != SWRAD_TS) {
			fail("The weathergenerator GWGEN requires: \n Tmax, Tmin, (Pressure & Specific Humidity) or rel.humidity, Windspeed, and SW radiation");
		}
		
		std::vector<double> mtemp;
		get_yearly_data(mtemp, spinup_temp, cf_temp, historic_timestep_temp);
		double xmtemp[12];
		for ( int i=0; i<12; i++) {
			xmtemp[i] = mtemp[i] - K2degC;
		}

		std::vector<double> mprec;
		get_yearly_data(mprec, spinup_prec, cf_prec, historic_timestep_prec);
		double xmprec[12];
		for ( int i=0; i<12; i++) {
			xmprec[i] = mprec[i];
		}

		std::vector<double> mwet;
		get_yearly_data(mwet, spinup_wetdays, cf_wetdays, historic_timestep_wetdays);
		double xmwet[12];
		for ( int i=0; i<12; i++) {
			xmwet[i] = mwet[i];
		}

		std::vector<double> minsol;
		get_yearly_data(minsol, spinup_insol, cf_insol, historic_timestep_insol);
		double xminsol[12];
		for ( int i=0; i<12; i++) {
			xminsol[i] = minsol[i];
		}
		std::copy(minsol.begin(),minsol.end(),xminsol);
		std::vector<double> mtmax;
		get_yearly_data(mtmax, spinup_max_temp, cf_max_temp, historic_timestep_max_temp);

		std::vector<double> mtmin;
		get_yearly_data(mtmin, spinup_min_temp, cf_min_temp, historic_timestep_min_temp);

		// Record shift of mean_temperature against mean of tmin/tmax for re-adjustment
		double xmdtr[12];
		double shift[12];
		for ( int i=0; i<12; i++) {
			xmdtr[i] = 0.5 * (mtmax[i] - mtmin[i]);
			shift[i] = 0.5 * (mtmax[i] + mtmin[i]) - mtemp[i];
		}

		std::vector<double> mwind;
		get_yearly_data(mwind, spinup_wind, cf_wind, historic_timestep_wind);
		double xmwind[12];
		for ( int i=0; i<12; i++) {
			xmwind[i] = mwind[i];
		}

		std::vector<double> mpres;
		if ( cf_pres ) {
			get_yearly_data(mpres, spinup_pres, cf_pres, historic_timestep_pres);
		}

		std::vector<double> mspecifichum;
		if ( cf_specifichum ) {
			get_yearly_data(mspecifichum, spinup_specifichum, cf_specifichum, historic_timestep_specifichum);
		}

		double xmrhum[12];
		std::vector<double> mrelhum;
		if ( cf_relhum ){ 
			get_yearly_data(mrelhum, spinup_relhum, cf_relhum, historic_timestep_relhum);
			for ( int i=0; i<12; i++) {
				xmrhum[i] = mrelhum[i];
			}
		}
		else if ( cf_pres && cf_specifichum ) {
			// compute rel. humidity if it can't be read from file
			for ( int i=0; i<12; i++) {
				xmrhum[i] = calc_relative_humidity(mtemp[i],mspecifichum[i],mpres[i]);
			}
		}

		// Use gwgen - correlated weather
		weathergen_get_met(gridcell,xmtemp,xmprec,xmwet,xminsol,xmdtr,
				   xmwind,xmrhum,dtemp,dprec,dinsol,ddtr,
				   dwind,drelhum);

		//produce tmin/tmax from daily temperature range plus shift
		int mon = 0;
		int accumday = 0;
		for (int i = 0; i < date.year_length(); ++i) {
			if ( i >= date.ndaymonth[mon]+accumday) {
				accumday += date.ndaymonth[mon];
				mon++;
			}
			dmin_temp[i] = dtemp[i] - 0.5 * ddtr[i];
			dmax_temp[i] = dtemp[i] + 0.5 * ddtr[i];
			// correct dmin and dmax against t_mean if available
			if ( cf_min_temp && cf_max_temp ) {
				dmin_temp[i] += shift[mon];
				dmax_temp[i] += shift[mon];
			}
		}
	}
	else {
		
		populate_daily_array(dtemp, spinup_temp, cf_temp, historic_timestep_temp, 0);
		populate_daily_prec_array(gridcell.seed);
		populate_daily_array(dinsol, spinup_insol, cf_insol, historic_timestep_insol, 0,
				     max_insolation(cf_standard_name_to_insoltype(cf_insol->get_standard_name())));
		
		if (cf_min_temp) {
			populate_daily_array(dmin_temp, spinup_min_temp, cf_min_temp, historic_timestep_min_temp, 0);
		}
		
		if (cf_max_temp) {
			populate_daily_array(dmax_temp, spinup_max_temp, cf_max_temp, historic_timestep_max_temp, 0);
		}
		
		if (cf_pres) {
			populate_daily_array(dpres, spinup_pres, cf_pres, historic_timestep_pres, 0);
		}
		
		if (cf_specifichum) {
			populate_daily_array(dspecifichum, spinup_specifichum, cf_specifichum, historic_timestep_specifichum, 0);
		}
		
		if (cf_wind) {
			populate_daily_array(dwind, spinup_wind, cf_wind, historic_timestep_wind, 0);
		}
		
		if (cf_relhum) {
			populate_daily_array(drelhum, spinup_relhum, cf_relhum, historic_timestep_relhum, 0);
		}
	}
	// Convert to units the model expects
	bool cloud_fraction_to_sunshine = (cf_standard_name_to_insoltype(cf_insol->get_standard_name()) == SUNSHINE);
	for (int i = 0; i < date.year_length(); ++i) {
		
		//is daily check for monthly outputs, avoid twice conversion
		if (is_daily(cf_temp) || weathergenerator == INTERP) {
			dtemp[i] -= K2degC;
			if (cf_min_temp) {
				dmin_temp[i] -= K2degC;
			}

			if (cf_max_temp) {
				dmax_temp[i] -= K2degC;
			}
		}
		
		if (cloud_fraction_to_sunshine) {
			// Invert from cloudiness to sunshine,
			// and convert fraction (0-1) to percent (0-100)
			dinsol[i] = (1-dinsol[i]) * 100.0;
		}
		
		//is daily check for monthly outputs, avoid twice conversion
		if (is_daily(cf_temp)) {
			if ((cf_pres && cf_specifichum) && !cf_relhum) {
				// compute relative humidity for BLAZE
				drelhum[i] = calc_relative_humidity(dtemp[i], dspecifichum[i], dpres[i]);
			}
		}
		else if ( firemodel == BLAZE && !cf_relhum ) {
			fail("BLAZE is switched on WITHOUT info on either specific humidity and pressure or relative humidity! \n" );
		}
	}

	// Move to next year in spinup dataset

	spinup_temp.nextyear();
	spinup_prec.nextyear();
	spinup_insol.nextyear();

	if (cf_wetdays) {
		spinup_wetdays.nextyear();
	}

	if (cf_min_temp) {
		spinup_min_temp.nextyear();
	}

	if (cf_max_temp) {
		spinup_max_temp.nextyear();
	}
	
	if (cf_pres) {
		spinup_pres.nextyear();
	}
	
	if (cf_specifichum) {
		spinup_specifichum.nextyear();
	}
	
	if (cf_wind) {
		spinup_wind.nextyear();
	}
	
	if (cf_relhum) {
		spinup_relhum.nextyear();
	}
	
	// Get monthly ndep values and convert to daily
	int ndep_year = date.get_calendar_year();
	// dprintf("Ndep_year = %d\n", ndep_year); //TA print year


	double mNHxdrydep[12], mNOydrydep[12];
	double mNHxwetdep[12], mNOywetdep[12];

	/*ndep.get_one_calendar_year(date.get_calendar_year(),
	                           mNHxdrydep, mNOydrydep,
							   mNHxwetdep, mNOywetdep);*/

	if (!(cf_mNHxdrydep && cf_mNOydrydep && cf_mNHxwetdep && cf_mNOywetdep)) {
		//dprintf("No N-deposition netCDF data used. Data of binary file is used.\n"); //TA print year
		ndep.get_one_calendar_year(date.get_calendar_year(),
			mNHxdrydep, mNOydrydep,
			mNHxwetdep, mNOywetdep);
	}
	else {
		int timestep = 0; // could just calculate the timestep (year-1850)*12
		GuessNC::CF::DateTime dd = cf_mNHxdrydep->get_date_time(timestep);
		GuessNC::CF::DateTime last_dd = cf_mNHxdrydep->get_date_time(cf_mNHxdrydep->get_timesteps() - 1);
		if (ndep_year > last_dd.get_year()) {
			ndep_year = last_dd.get_year();
		}
		if (ndep_year > dd.get_year()) { // advance timestep to Jan in the year
			while (dd.get_year() < ndep_year) {
				timestep++;
				dd = cf_mNHxdrydep->get_date_time(timestep);
			}
		}

		// read 12 monthly values
		for (int m = 0; m < 12; m++) {
			//dprintf("Ndep_year = %d\n", ndep_year); //TA print year
			mNHxdrydep[m] = cf_mNHxdrydep->get_value(timestep); // NO convert kg m-2 s-1 -> kg m-2 d-1
			mNOydrydep[m] = cf_mNOydrydep->get_value(timestep); //
			mNHxwetdep[m] = cf_mNHxwetdep->get_value(timestep); //
			mNOywetdep[m] = cf_mNOywetdep->get_value(timestep); //
			timestep++;
			//dprintf("mNHxdrydep = %f\n", mNHxdrydep[1]); //TA print mNHxdrydep of first month
		}
	}


	// Divide pdep into dry and wet
	double mpdrydep[12], mpwetdep[12];

	for (int m = 0; m < 12; m++) {
		mpdrydep[m] = 0.0;
		mpwetdep[m] = 0.0;
	}

	// Phosphorus deposition
	if (!(cf_mPOxdrydep && cf_mPOxwetdep)) {
		/*double gridcell_mpdep[12];
		get_monthly_pdep(gridcell.get_lat(), gridcell.get_lon(), gridcell_mpdep);
		for (int m = 0; m < 12; m++) {
			mpdrydep[m] = gridcell_mpdep[m] / 2.0;
			mpwetdep[m] = gridcell_mpdep[m] / 2.0;
		}*/
		if (date.day == 0 && date.year == 0) {
			get_monthly_pdep(gridcell.get_lat(), gridcell.get_lon(), gridcell.climate.mpdep);
			for (int m = 0; m < 12; m++) {
				mpdrydep[m] = gridcell.climate.mpdep[m] / 2.0;
				mpwetdep[m] = gridcell.climate.mpdep[m] / 2.0;
			}
		}
	}
	else {
		int timestep = 0; // could just calculate the timestep (year-1850)*12
		GuessNC::CF::DateTime dd = cf_mPOxdrydep->get_date_time(timestep);
		GuessNC::CF::DateTime last_dd = cf_mPOxdrydep->get_date_time(cf_mPOxdrydep->get_timesteps() - 1);
		if (ndep_year > last_dd.get_year()) {
			ndep_year = last_dd.get_year();
		}
		if (ndep_year > dd.get_year()) { // advance timestep to Jan in the year
			while (dd.get_year() < ndep_year) {
				timestep++;
				dd = cf_mPOxdrydep->get_date_time(timestep);
			}
		}

		// read 12 monthly values
		for (int m = 0; m < 12; m++) {
			//dprintf("Ndep_year = %d\n", ndep_year); //TA print year
			mpdrydep[m] = cf_mPOxdrydep->get_value(timestep); // NO convert kg m-2 s-1 -> kg m-2 d-1
			mpwetdep[m] = cf_mPOxwetdep->get_value(timestep); //
			timestep++;
			//dprintf("mNHxdrydep = %f\n", mNHxdrydep[1]); //TA print mNHxdrydep of first month
		}
	}

	// Distribute N deposition
	distribute_ndep(mNHxdrydep, mNOydrydep,
					mNHxwetdep, mNOywetdep,
					dprec,dNH4dep,dNO3dep);

	// Distribute P deposition
	distribute_pdep(mpdrydep, mpwetdep,
		dprec, dpdep);

	// Phosphorus weathering
	if (param["file_pwtr"].str != "") {
		if (date.day == 0 && date.year == 0)
			get_yearly_pwtr_params(gridcell.get_lat(), gridcell.get_lon(), gridcell.climate.pwtr_bi,
				gridcell.climate.pwtr_pcont, gridcell.climate.pwtr_shield, gridcell.climate.pwtr_ea);
	}
}

void CFInput::getlandcover(Gridcell& gridcell) {

	landcover_input.getlandcover(gridcell);
	landcover_input.get_land_transitions(gridcell);
}

bool CFInput::getclimate(Gridcell& gridcell) {

	Climate& climate = gridcell.climate;

	GuessNC::CF::DateTime last_date = last_day_to_simulate(cf_temp);

	if (later_day(date, last_date)) {
		++current_gridcell;
		return false;
	}

	climate.co2 = co2[date.get_calendar_year()];

	if (date.day == 0) {
		populate_daily_arrays(gridcell);
	}

	climate.temp   = dtemp[date.day];
	climate.prec   = dprec[date.day];
	climate.insol  = dinsol[date.day];
	climate.relhum = drelhum[date.day];
	climate.u10    = dwind[date.day];
	climate.tmax   = dmax_temp[date.day];
	climate.tmin   = dmin_temp[date.day];
	climate.dtr    = ddtr[date.day];

	// Nitrogen deposition
	gridcell.dNH4dep = dNH4dep[date.day];
	gridcell.dNO3dep = dNO3dep[date.day];

	// Phosphorus deposition
	gridcell.dpdep = dpdep[date.day];

	// bvoc
	if(ifbvoc){
		if (cf_min_temp && cf_max_temp) {
			climate.dtr = dmax_temp[date.day] - dmin_temp[date.day];
		}
		else {
			fail("When BVOC is switched on, valid paths for minimum and maximum temperature must be given.");
		}
	}

	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {

			int first_historic_year = cf_temp->get_date_time(0).get_year();
			int last_historic_year = cf_temp->get_date_time(cf_temp->get_timesteps()-1).get_year();
			int historic_years = last_historic_year - first_historic_year + 1;

			int years_to_simulate = nyear_spinup + historic_years;

			int cells_done = (int)distance(gridlist.begin(), current_gridcell);

			double progress = (double)(cells_done*years_to_simulate+date.year)/
				(gridlist.size()*(double)years_to_simulate);
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;

}

void CFInput::load_spinup_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
                               GenericSpinupData& spinup_data) {

	const std::string error_message =
		format_string("Not enough data to build spinup, at least %d years needed",
		              NYEAR_SPINUP_DATA);

	GenericSpinupData::RawData source;

	int timestep = 0;

	bool daily = is_daily(cf_var);
	// for now, assume that each data set is either daily or monthly
	bool monthly = !daily;

	// Skip the first year if data doesn't start at the beginning of the year
	while ((daily && !first_day_of_year(cf_var->get_date_time(timestep))) ||
	       (monthly && !first_month_of_year(cf_var->get_date_time(timestep)))) {
		++timestep;

		if (timestep >= cf_var->get_timesteps()) {
			fail(error_message.c_str());
		}
	}

	// Get all the values for the first NYEAR_SPINUP_DATA years,
	// and put them into source
	for (int i = 0; i < NYEAR_SPINUP_DATA; ++i) {
		std::vector<double> year(daily ? GenericSpinupData::DAYS_PER_YEAR : 12);

		for (size_t i = 0; i < year.size(); ++i) {
			if (timestep < cf_var->get_timesteps()) {
				GuessNC::CF::DateTime dt = cf_var->get_date_time(timestep);

				if (daily && dt.get_month() == 2 && dt.get_day() == 29) {
					++timestep;
				}
			}

			if (timestep >= cf_var->get_timesteps()) {
				fail(error_message.c_str());
			}

			year[i] = cf_var->get_value(timestep);
			++timestep;
		}

		source.push_back(year);
	}

	spinup_data.get_data_from(source);
}

namespace {

// Help function for call to remove_if below - checks if a pointer is null
bool is_null(const GuessNC::CF::GridcellOrderedVariable* ptr) {
	return ptr == 0;
}

}

std::vector<GuessNC::CF::GridcellOrderedVariable*> CFInput::all_variables() const {
	std::vector<GuessNC::CF::GridcellOrderedVariable*> result;
	result.push_back(cf_temp);
	result.push_back(cf_prec);
	result.push_back(cf_insol);
	result.push_back(cf_wetdays);
	result.push_back(cf_min_temp);
	result.push_back(cf_max_temp);
	result.push_back(cf_pres);
	result.push_back(cf_specifichum);
	result.push_back(cf_wind);
	result.push_back(cf_relhum);

	// Get rid of null pointers
	result.erase(std::remove_if(result.begin(), result.end(), is_null),
	             result.end());

	return result;
}

#endif // HAVE_NETCDF
