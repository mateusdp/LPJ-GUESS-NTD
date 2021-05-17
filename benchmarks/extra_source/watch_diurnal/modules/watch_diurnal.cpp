///////////////////////////////////////////////////////////////////////////////////////
/// \file watch_diurnal.cpp
/// \brief Input module for reading in WATCH diurnal data from NetCDF files
///
/// \TODO This input module have not implemented reading in of relative humidity and 
/// windspeed (needed) for Blaze, currently taken from cru_input. This has the implication 
/// that the climate variables are not synced. 
///
/// $Date: 2019-10-28 18:48:52 +0100 (Mon, 28 Oct 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "watch_diurnal.h"
#include <netcdf.h>
#include <assert.h>
#include "guess.h"
#include "driver.h"

REGISTER_INPUT_MODULE("watch_diurnal", WATCHDiurnalInput)

void handle_error(int status, const char* message) {
	if (status != NC_NOERR) {
		fail("NetCDF error: %s: %s\n", nc_strerror(status), message);
	}
}

int open_ncdf(const char* fname) {
	int netcdf_id;
	int status = nc_open(fname, NC_NOWRITE, &netcdf_id);
	handle_error(status, (std::string("Cannot open NetCDF file: ") + fname).c_str());
	return netcdf_id;
}

void load_gridlist(const char* dir_name, std::vector<landpoint>& landpoints) {
	using std::vector;
	using std::string;

	string gridlist_file = string(dir_name) + "/WFD-land-lat-long-z.nc";
	int ncid = open_ncdf(gridlist_file.c_str());

	// Get the land dimension and figure out how many grid cells there are

	int landid;
	int status = nc_inq_dimid(ncid, "land", &landid);
	handle_error(status, "land dimension");

	size_t num_gridcells;
	status = nc_inq_dimlen(ncid, landid, &num_gridcells);
	handle_error(status, "land dimension length");

	// Read in all longitudes and latitudes into two vectors

	vector<double> lons(num_gridcells), lats(num_gridcells);

	int lonid, latid;
	status = nc_inq_varid(ncid, "Longitude", &lonid);
	handle_error(status, "Longitude variable id");

	status = nc_inq_varid(ncid, "Latitude", &latid);
	handle_error(status, "Latitude variable id");

	nc_get_var_double(ncid, lonid, &lons.front());
	nc_get_var_double(ncid, latid, &lats.front());

	// Transfer the coordinates to the landpoints vector

	landpoints.resize(num_gridcells);

	for (size_t i = 0; i < landpoints.size(); i++) {
		landpoints[i].first = lons[i];
		landpoints[i].second = lats[i];
	}

	nc_close(ncid);
}

int get_cell_id(double lon, double lat, const std::vector<landpoint>& landpoints) {
	for (size_t i = 0; i < landpoints.size(); i++) {
		if (landpoints[i].first == lon && landpoints[i].second == lat) {
			return i;
		}
	}
	fail("Cell (%g,%g) wasn't found in the NetCDF grid\n", lon, lat);
	return -1;
}

void load_watch_data(const char* dir_name,
                     const char* var_name,
                     int cell_id,
                     std::vector<double>& data,
                     bool diurnal) {
	using std::vector;

	// Open the NetCDF file for this variable and grid cell

	xtring fname;
	fname.printf("%s/%s/gc_%i.nc", dir_name, var_name, cell_id);

	int ncid = open_ncdf(fname);

	// Figure out the number of time steps

	int tstep_id;
	int status = nc_inq_dimid(ncid, "tstep", &tstep_id);
	handle_error(status, "tstep dimension");

	size_t num_timesteps;
	status = nc_inq_dimlen(ncid, tstep_id, &num_timesteps);
	handle_error(status, "tstep dimension length");

	// Read in all the data from the variable

	int var_id;
	status = nc_inq_varid(ncid, var_name, &var_id);
	handle_error(status, var_name);

	vector<double> raw_data(num_timesteps);
	status = nc_get_var_double(ncid, var_id, &raw_data.front());
	handle_error(status, var_name);

	// Transfer data from 'raw_data' to 'data', skipping leap days

	data.clear();
	data.reserve(num_timesteps);

	vector<double>::const_iterator raw_data_itr = raw_data.begin();

	const size_t steps_per_day = diurnal ? SUBDAILY : 1;

	int calendar_year = FIRST_WATCH_YEAR;

	// loop through the whole raw_data vector
	while (raw_data_itr != raw_data.end()) {

		// copy one year in each iteration
		for (int i = 0; i < 365; ++i) {
			if (i == 59 && Date::is_leap(calendar_year)) {
				// leap day, skip it
				raw_data_itr += steps_per_day;
			}
			data.insert(data.end(), raw_data_itr, raw_data_itr + steps_per_day);
			raw_data_itr += steps_per_day;
		}
		++calendar_year;
	}

	// Make sure we agree on leap days etc.
	assert(raw_data_itr == raw_data.end());

	nc_close(ncid);
}


WATCHDiurnalInput::WATCHDiurnalInput()
	: diurnal(false) {

	declare_parameter("diurnal", &diurnal, "If specified, diurnal version will be run (0,1)");
}

void WATCHDiurnalInput::init() {
	CRUInput::init();

	watch_dir = param["watch_dir"].str;
	load_gridlist(watch_dir, landpoints);
}

bool WATCHDiurnalInput::getgridcell(Gridcell& gridcell) {
	if (CRUInput::getgridcell(gridcell)) {
		// Load WATCH subdaily variables
		int cell_id = get_cell_id(gridcell.get_lon(), gridcell.get_lat(), landpoints);

		load_watch_data(watch_dir, "Tair",   cell_id, watch_temp,   true);
		load_watch_data(watch_dir, "SWdown", cell_id, watch_swdown, true);
		load_watch_data(watch_dir, "Rainf",  cell_id, watch_rainf,  false);
		load_watch_data(watch_dir, "Snowf",  cell_id, watch_snowf,  false);

		// Get nitrogen deposition data.
		/* Since the historic data set does not reach decade 2010-2019,
		 * we need to use the RCP data for the last decade. */
		ndep.getndep(param["file_ndep"].str, gridcell.get_lon(), gridcell.get_lat(), Lamarque::RCP60);

		// The insolation data will be sent (in function getclimate, below)
		// as incoming shortwave radiation, averages are over 24 hours

		gridcell.climate.instype = SWRAD_TS;

		return true;
	}
	else {
		return false;
	}
}

bool WATCHDiurnalInput::getclimate(Gridcell& gridcell) {
	if (CRUInput::getclimate(gridcell)) {

		Climate& climate = gridcell.climate;

		if (date.day == 0) {
			double dprec[365];
	        double mNHxdrydep[12], mNOydrydep[12], mNHxwetdep[12], mNOywetdep[12];
			ndep.get_one_calendar_year(date.year - nyear_spinup + FIRST_WATCH_YEAR,
			                           mNHxdrydep, mNOydrydep,
									   mNHxwetdep, mNOywetdep);

			//get_monthly_ndep(FIRST_WATCH_YEAR + date.year - nyear_spinup, mndrydep, mnwetdep);

			// Distribute N deposition - without rain days
			std::fill_n(dprec, 365, 0);
			// Distribute N deposition
			distribute_ndep(mNHxdrydep, mNOydrydep,
							mNHxwetdep, mNOywetdep,
							dprec,dNH4dep,dNO3dep);
		}

		int year = date.year < nyear_spinup ?
			date.year % NYEAR_SPINUP_DATA : date.year - nyear_spinup;

		size_t daily_index = year * 365 + date.day;
		size_t subdaily_index = daily_index * SUBDAILY;
		size_t subdaily_end = subdaily_index + SUBDAILY;

		if (daily_index >= watch_rainf.size()) {
			// no more forcing data left, so we're finished with this grid cell
			return false;
		}

		climate.temps.assign(watch_temp.begin()+subdaily_index,
		                     watch_temp.begin()+subdaily_end);
		climate.insols.assign(watch_swdown.begin()+subdaily_index,
		                      watch_swdown.begin()+subdaily_end);

		// Convert temperatures (K -> C)
		for (size_t i = 0; i < SUBDAILY; ++i) {
			climate.temps[i] -= K2degC;
		}

		climate.temp = mean(&climate.temps.front(), SUBDAILY);
		climate.insol = mean(&climate.insols.front(), SUBDAILY);
		climate.prec = (watch_rainf[daily_index] + watch_snowf[daily_index]);
		climate.prec *= 24 * 3600; // mm/s -> mm/day

		if (ifbvoc && !diurnal) {
			std::vector<double>::const_iterator start, end;
			start = climate.temps.begin();
			end = climate.temps.end();
			climate.dtr = std::max_element(start, end) - std::min_element(start, end);
		}
		date.subdaily = diurnal ? SUBDAILY : 1;

		// Nitrogen deposition
		gridcell.dNH4dep = dNH4dep[date.day];
		gridcell.dNO3dep = dNO3dep[date.day];

		return true;
	}
	else {
		return false;
	}
}
