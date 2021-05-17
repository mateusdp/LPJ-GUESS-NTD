///////////////////////////////////////////////////////////////////////////////////////
/// \file watch_diurnal.h
/// \brief Input module for reading in WATCH diurnal data from NetCDF files
///
/// $Date: 2019-10-28 18:48:52 +0100 (Mon, 28 Oct 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_WATCH_DIURNAL_H
#define LPJ_GUESS_WATCH_DIURNAL_H

#ifndef HAVE_NETCDF
#error "NetCDF library not found by build system!"
#endif

#include "cruinput.h"

const int SUBDAILY = 8;
const int FIRST_WATCH_YEAR = 1901;

/// A lon,lat pair, used by load_gridlist function below
typedef std::pair<double, double> landpoint;

/// Help function to deal with status codes from NetCDF library
void handle_error(int status, const char* message);

/// Opens a NetCDF file and returns its id
int open_ncdf(const char* fname);

/// Reads in all the WATCH coordinates into a vector.
/** The index of a coordinate in the landpoints vector is the grid cell's id,
 *  which is used to find the NetCDF file for the grid cell.
 */
void load_gridlist(const char* dir_name, std::vector<landpoint>& landpoints);

/// Returns the grid cell id for a coordinate
int get_cell_id(double lon, double lat, const std::vector<landpoint>& landpoints);

/// Read in all data for a single variable
/** Data for leap years is skipped.
 *
 * \param dir_name Path to WATCH directory with all NetCDF files
 * \param var_name Which variable to read in
 * \param cell_id  id number for the grid cell to read in
 * \param diurnal  Whether the variable has subdaily or daily data */
void load_watch_data(const char* dir_name,
                     const char* var_name,
                     int cell_id, 
                     std::vector<double>& data,
                     bool diurnal);

/// Input module for reading WATCH diurnal data from NetCDF files
class WATCHDiurnalInput : public CRUInput {
public:

	WATCHDiurnalInput();

	// Reimplemented functions from CRUInput

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

private:
	/// Whether to run in diurnal mode or not
	bool diurnal;

	/// Directory of the WATCH NetCDF files
	xtring watch_dir;
	
	/// Used to find grid cell ids, given a coordinate
	std::vector<landpoint> landpoints;
	
	/// WATCH forcing data
	/** Flat arrays taken straight from the NetCDF files, same units,
	 *  but leap days taken out. Read in for the current grid cell
	 *  in getgridcell().
	 */
	std::vector<double> watch_temp, watch_swdown, watch_rainf, watch_snowf;

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;

	/// Daily N deposition for current year
	double dNH4dep[Date::MAX_YEAR_LENGTH];
	double dNO3dep[Date::MAX_YEAR_LENGTH];
};

#endif // LPJ_GUESS_WATCH_DIURNAL_H
