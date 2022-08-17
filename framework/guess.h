///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.h
/// \brief Framework header file, LPJ-GUESS Combined Modular Framework
///
/// This header file contains:
///  (1) definitions of all main classes used by the framework and modules. Modules may
///      require classes to contain certain member variables and functions (see module
///      source files for details).
///  (2) other type, constant and function definitions to be accessible throughout the
///      model code.
///  (3) a forward declaration of the framework function if this is not the main
///      function.
///
/// \author Ben Smith
/// $Date: 2022-09-13 10:47:57 +0200 (Tue, 13 Sep 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_GUESS_H
#define LPJ_GUESS_GUESS_H

///////////////////////////////////////////////////////////////////////////////////////
// #INCLUDES FOR LIBRARY HEADER FILES
// C/C++ libraries required for member functions of classes defined in this file.
// These libraries will also be available globally (so omit these #includes from source
// files). In addition to various standard C/C++ runtime libraries, the framework
// requires the following libraries (individual modules may use additional libraries)
//
// GUTIL
//   Includes class xtring, providing functionality for pointer-free dynamic handling
//   of character strings; wherever possible in LPJ-GUESS, strings are represented as
//   objects of type xtring rather than simple arrays of type char. GUTIL also provides
//   templates for dynamic collection classes (list arrays of various types), argument
//   processing for printf-style functions, timing functions and other utilities.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "gutil.h"
#include <vector>
#include <algorithm>
#include "shell.h"
#include "guessmath.h"
#include "archive.h"
#include "parameters.h"
#include "guesscontainer.h"
#include "soil.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL ENUMERATED TYPE DEFINITIONS

/// Life form class for PFTs (trees, grasses)
typedef enum {NOLIFEFORM, TREE, GRASS, MOSS} lifeformtype;

/// Phenology class for PFTs
typedef enum {NOPHENOLOGY, EVERGREEN, RAINGREEN, SUMMERGREEN, CROPGREEN, ANY} phenologytype;

/// Biochemical pathway for photosynthesis (C3 or C4)
typedef enum {NOPATHWAY, C3, C4} pathwaytype;

/// Leaf physiognomy types for PFTs
typedef enum {NOLEAFTYPE, NEEDLELEAF, BROADLEAF} leafphysiognomytype;

/// The level of verbosity of LPJ-GUESS. Decides the amount of information that is written to the log-file.
typedef enum {ERROR, WARNING, INFO, DEBUG_WARNING} verbositylevel;

/// Units for insolation driving data
/** Insolation can be expressed as:
 *
 *  - Percentage sunshine
 *  - Net instantaneous downward shortwave radiation flux (W/m2)
 *  - Total (i.e. with no correction for surface albedo) instantaneous downward
 *    shortwave radiation flux (W/m2)
 *
 *  Radiation flux can be interpreted as W/m2 during daylight hours, or averaged
 *  over the whole time step which it represents (24 hours in daily mode). For
 *  this reason there are two enumerators for these insolation types (e.g. SWRAD
 *  and SWRAD_TS).
 */
typedef enum {
	/// No insolation type chosen
	NOINSOL,
	/// Percentage sunshine
	SUNSHINE,
	/// Net shortwave radiation flux during daylight hours (W/m2)
	NETSWRAD,
	/// Total shortwave radiation flux during daylight hours (W/m2)
	SWRAD,
	/// Net shortwave radiation flux during whole time step (W/m2)
	NETSWRAD_TS,
	/// Total shortwave radiation flux during whole time step (W/m2)
	SWRAD_TS
} insoltype;

/// CENTURY pool names, NSOMPOOL number of SOM pools
typedef enum {SURFSTRUCT, SOILSTRUCT, SOILMICRO, SURFHUMUS, SURFMICRO, SURFMETA, SURFFWD, SURFCWD,
	SOILMETA, SLOWSOM, PASSIVESOM, NSOMPOOL} pooltype;

/// Irrigation type for PFTs
typedef enum {RAINFED, IRRIGATED} hydrologytype;
/// Intercrop type for PFTs
typedef enum {NOINTERCROP, NATURALGRASS} intercroptype;

/// Seasonality type of gridcell
/** 0:SEASONALITY_NO			No seasonality
 *  1:SEASONALITY_PREC			Precipitation seasonality only
 *  2:SEASONALITY_PRECTEMP		Both temperature and precipitation seasonality, but "weak" temperature seasonality (coldest month > 10degC)
 *  3:SEASONALITY_TEMP			Temperature seasonality only
 *  4:SEASONALITY_TEMPPREC		Both temperature and precipitation seasonality, but temperature most important (coldest month < 10degC)
 *  5:SEASONALITY_TEMPWARM		Temperature seasonality, always above 10 degrees (currently not used)
 */
typedef enum {SEASONALITY_NO, SEASONALITY_PREC, SEASONALITY_PRECTEMP, SEASONALITY_TEMP, SEASONALITY_TEMPPREC} seasonality_type;

/// Precipitation seasonality type of gridcell
/** 0:DRY						(minprec_pet20<=0.5 && maxprec_pet20<=0.5)
 *  1:DRY_INTERMEDIATE			(minprec_pet20<=0.5 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
 *  2:DRY_WET					(minprec_pet20<=0.5 && maxprec_pet20>1.0)
 *  3:INTERMEDIATE				(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
 *  4:INTERMEDIATE_WET			(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>1.0)
 *  5:WET						(minprec_pet20>1.0 && maxprec_pet20>1.0)
 */
typedef enum {DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET} prec_seasonality_type;


/// Temperature seasonality type of gridcell
/** 0:COLD						(mtemp_max20<=10)
 *  1:COLD_WARM					(mtemp_min20<=10 && mtemp_max20>10 && mtemp_max20<=30)
 *  2:COLD_HOT					(mtemp_min20<=10 && mtemp_max20>30)
 *  3:WARM						(mtemp_min20>10 && mtemp_max20<=30)
 *  4:WARM_HOT					(mtemp_min20>10 && mtemp_max20>30)
 *  5:HOT						(mtemp_min20>30)
 */
typedef enum {COLD, COLD_WARM, COLD_HOT, WARM, WARM_HOT, HOT} temp_seasonality_type;

/// Gas type (used in methane code)
/** 
  */
typedef enum {O2gas, CO2gas, CH4gas} gastype;

/// Nitrogen preferance
typedef enum {NO, NH4, NO3} n_pref_type;
///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS

/// number  of soil layers modelled
const int NSOILLAYER_UPPER = 5;
const int NSOILLAYER_LOWER = NSOILLAYER - NSOILLAYER_UPPER;

/// bvoc: number of monoterpene species used
const int NMTCOMPOUNDS=NMTCOMPOUNDTYPES;

// SOIL DEPTH VALUES

/// soil upper layer depth (mm)
const double SOILDEPTH_UPPER = 500.0;
/// soil lower layer depth (mm)
const double SOILDEPTH_LOWER = 1000.0;

/// Depth of sublayer at top of upper soil layer, from which evaporation is
//  possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
//  Must be a multiple of Dz_soil
const double SOILDEPTH_EVAP = 200.0;

/// Year at which to calculate equilibrium soil carbon
const int SOLVESOM_END=400;

/// Year at which to begin documenting means for calculation of equilibrium soil carbon
const int SOLVESOM_BEGIN = 350;

/// Number of years to average growth efficiency over in function mortality
const int NYEARGREFF = 5;

/// Coldest day in N hemisphere (January 15)
/** Used to decide when to start counting GDD's and leaf-on days
 *  for summergreen phenology.
 */
const int COLDEST_DAY_NHEMISPHERE = 14;

/// Coldest day in S hemisphere (July 15)
/** Used to decide when to start counting GDD's and leaf-on days
 *  for summergreen phenology.
 */
const int COLDEST_DAY_SHEMISPHERE = 195;

/// Warmest day in N hemisphere (same as COLDEST_DAY_SHEMISPHERE)
const int WARMEST_DAY_NHEMISPHERE = COLDEST_DAY_SHEMISPHERE;

/// Warmest day in S hemisphere (same as COLDEST_DAY_NHEMISPHERE)
const int WARMEST_DAY_SHEMISPHERE = COLDEST_DAY_NHEMISPHERE;

/// number of years to average aaet over in function soilnadd
const int NYEARAAET = 5;

/// number of years to average max snow depth over in function soilnadd
const int NYEARMAXSNOW = 20;

/// Priestley-Taylor coefficient (conversion factor from equilibrium evapotranspiration to PET)
const double PRIESTLEY_TAYLOR = 1.32;

// Solving Century SOM pools

/// fraction of nyear_spinup minus freenyears at which to begin documenting for calculation of Century equilibrium
const double SOLVESOMCENT_SPINBEGIN  = 0.1;
/// fraction of nyear_spinup minus freenyears at which to end documentation and start calculation of Century equilibrium
const double SOLVESOMCENT_SPINEND    = 0.3;

/// Kelvin to deg C conversion
const double K2degC = 273.15;

/// Maximum number of crop rotation items
const int NROTATIONPERIODS_MAX = 3;

/// Conversion factor for CO2 from ppmv to mole fraction
const double CO2_CONV = 1.0e-6;

/// Initial carbon allocated to crop organs at sowing, kg m-2
const double CMASS_SEED = 0.01;

/// Precision in land cover fraction input
const double INPUT_PRECISION = 1.0e-14;
const double INPUT_ERROR = 0.5e-6;
const double INPUT_RESOLUTION = INPUT_PRECISION - INPUT_PRECISION * INPUT_ERROR;

/// Averaging interval for average maximum annual fapar (SIMFIRE)
const int AVG_INTERVAL_FAPAR = 3;

/// Averaging interval for biome averaging (SIMFIRE)
const int N_YEAR_BIOMEAVG = 3;

///////////////////////////////////////////////////////////////////////////////////////
// FORWARD DECLARATIONS OF CLASSES DEFINED IN THIS FILE
// Forward declarations of classes used as types (e.g. for reference variables in some
// classes) before they are actually defined

class Date;
class Stand;
class Patch;
class Vegetation;
class Individual;
class Gridcell;
class Patchpft;

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are defined in the framework source code file, and are accessible
// throughout the code

/// Object describing timing stage of simulation
extern Date date;

/// Number of possible PFTs
extern int npft;
/// Number of stand types in stlist
extern int nst;
/// Number of stand types per land cover
extern int nst_lc[NLANDCOVERTYPES];
/// Number of management types in stlist
extern int nmt;

/// Routine for handling slightly out of bounds water contents arising in hydrology and soil water thaw/freezing calculations.
inline void oob_check_wcont(double &wc_in) {

	// Minimum nonzero wcont allowed
	// This function will remove very small water amounts from wcont variables and adjust to 1 if we exceed it by tiny amounts
	const double min_mm = 0.000000001;

	if (DEBUG_SOIL_WATER) {
		if (wc_in <= -1.0 * min_mm || wc_in >= 1.0 + min_mm) {
			fail("Out of bounds detected in oob_check_wcont - water content arising in hydrology and soil water thaw/freezing calculations! wcont %g\n", wc_in);
		}
	}

	if (wc_in != 0.0 && wc_in < min_mm && wc_in > -1.0 * min_mm) {
		wc_in = 0.0;
	}

	if (wc_in > 1.0 && wc_in < 1.0 + min_mm) {
		wc_in = 1.0;
	}
}

/// General purpose object for handling simulation timing.
/** In general, frameworks should use a single Date object for all simulation
 *  timing.
 *
 *  Member variables of the class (see below) provide various kinds of calender
 *  and timing information, assuming init has been called to initialise the
 *  object, and next() has been called at the end of each simulation day.
 */
class Date {

	// MEMBER VARIABLES

public:

	/// Maximum number of days in an LPJ-GUESS simulation year
	/** The standard version doesn't yet support leap years. */
	static const int MAX_YEAR_LENGTH = 365;

	/// number of days in each month (0=January - 11=December)
	int ndaymonth[12];

	/// julian day of year (0-364; 0=Jan 1)
	int day;

	/// day of current month (0=first day)
	int dayofmonth;

	/// month number (0=January - 11=December)
	int month;

	/// year since start of simulation (0=first simulation year)
	int year;

	/// number of subdaily periods in a day (to be set in IO module)
	int subdaily;

	/// julian day for middle day of each month
	int middaymonth[12];

	/// true if last year of simulation, false otherwise
	bool islastyear;

	/// true if last month of year, false otherwise
	bool islastmonth;

	/// true if last day of month, false otherwise
	bool islastday;

	/// true if middle day of month, false otherwise
	bool ismidday;

	/// The calendar year corresponding to simulation year 0
	int first_calendar_year;

private:

	int nyear;

	// MEMBER FUNCTIONS

public:

	/// Constructor function called automatically when Date object is created
	/** Do not call explicitly. Initialises some member variables. */
	Date() {
		const int data[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
		int month;
		int dayct = 0;
		for (month=0; month<12; month++) {
			ndaymonth[month] = data[month];
			middaymonth[month] = dayct + data[month] / 2;
			dayct += data[month];
		}
		subdaily = 1;
		first_calendar_year = 0;
	}

	/// Initialises date to day 0 of year 0 and sets intended number of simulation years
	/** Intended number of simulation years is only used to set islastyear flag,
	 *  actual simulation may be longer or shorter.
	 *
	 *  \param nyearsim  Intended number of simulation years
	 */
	void init(int nyearsim)	{
		nyear = nyearsim;
		day = month=year = dayofmonth = 0;
		islastmonth = islastday = ismidday = false;
		if (nyear == 1) islastyear = true;
		else islastyear = false;
	}

	/// Call at end of every simulation day to update member variables.
	void next() {
		if (islastday) {
			if (islastmonth) {
				dayofmonth = 0;
				day = 0;
				month = 0;
				year++;
				if (year == nyear - 1) islastyear = true;
				islastmonth = false;
			}
			else {
				day++;
				dayofmonth = 0;
				month++;
				if (month == 11) islastmonth = true;
			}
			islastday = false;
		}
		else {
			day++;
			dayofmonth++;
			if (dayofmonth == ndaymonth[month] / 2) ismidday = true;
			else {
				ismidday = false;
				if (dayofmonth == ndaymonth[month] - 1) islastday = true;
			}
		}
	}

	// \returns index (0-11) of previous month (11 if currently month 0).
	int prevmonth() {
		if (month > 0) return month - 1;
		return 11;
	}

	/// \returns index of next month (0 if currently month 11)
	int nextmonth() {
		if (month < 11) return month+1;
		return 0;
	}

	/// Check if the year is leap
	/** \param year        Calendar year
	*   The algorith is as follows: only year that are divisible by 4 could
	*   potentially be leap (e.g., 1904), however, not if they're divisble by
	*   100 (e.g., 1900 is not leap), unless they're divisble by 400 (e.g., 2000
	*   is still leap).
	*/
	static bool is_leap(int year) {
		return (!(year % 4) && (year % 100 | !(year % 400)));
	}

	/// Whether the current mode is diurnal
	bool diurnal() const { return subdaily > 1; }

	/// Sets calendar year for simulation year 0
	/** Astronomical year numbering is used, so year 1 BC is represented by 0,
	 *  2 BC = -1 etc. See ISO 8601.
	 */
	void set_first_calendar_year(int calendar_year) {
		first_calendar_year = calendar_year;
	}

	/// Returns the calendar year corresponding to the current simulation year
	/** Astronomical year numbering is used, so year 1 BC is represented by 0,
	 *  2 BC = -1 etc. See ISO 8601.
	 */
	int get_calendar_year() const {
		return year + first_calendar_year;
	}

	/// \returns the number of days in the current simulation year
	/** For this function to work properly in simulations with varying number
	 *  of days per year, the set_first_calendar_year must have been called first.
	 *
	 *  Currently there is no support for leap years, so this function
	 *  always returns 365. */
	int year_length() const {
		return MAX_YEAR_LENGTH;
	}

	/// Step n days from a date.
	/** Current implementation does not consider leap days, and the same 
	  * apply for the current use of the function through-out the model.
      */
	static int stepfromdate(int day, int step) {

		if(day < 0)			// a negative value should not be a valid day
			return -1;
		else if(day + step > 0)
			return (day + step) % MAX_YEAR_LENGTH;
		else if(day + step < 0)
			return day + step + MAX_YEAR_LENGTH;
		else
			return 0;
	}
};

/// Object describing sub-daily periods
class Day {
public:
	/// Whether sub-daily period first/last within day (both true in daily mode)
	bool isstart, isend;

	/// Ordinal number of the sub-daily period [0, date.subdaily)
	int period;

	/// Constructs beginning of the day period (the only one in daily mode)
	Day() {
		isstart = true;
		isend = !date.diurnal();
		period = 0;
	}

	/// Advances to the next sub-daily period
	void next() {
		period++;
		isstart = false;
		isend = period == date.subdaily - 1;
	}
};

/// Object updating gridcell mass balance; currently used in framework()
class MassBalance : public Serializable  {

	int start_year;
	double ccont;
	double ccont_zero;
	double ccont_zero_scaled;
	double cflux;
	double cflux_zero;

	double ncont;
	double ncont_zero;
	double ncont_zero_scaled;
	double nflux;
	double nflux_zero;

	double pcont;
	double pcont_zero;
	double pcont_zero_scaled;
	double pflux;
	double pflux_zero;

public:
	MassBalance() {

		start_year = nyear_spinup;
		ccont = 0.0;
		ccont_zero = 0.0;
		ccont_zero_scaled = 0.0;
		cflux = 0.0;
		cflux_zero = 0.0;
		ncont = 0.0;
		ncont_zero = 0.0;
		ncont_zero_scaled = 0.0;
		nflux = 0.0;
		nflux_zero = 0.0;
		pcont = 0.0;
		pcont_zero = 0.0;
		pcont_zero_scaled = 0.0;
		pflux = 0.0;
		pflux_zero = 0.0;
	}

	MassBalance(int start_yearX) {

		start_year = start_yearX;
		ccont = 0.0;
		ccont_zero = 0.0;
		ccont_zero_scaled = 0.0;
		cflux = 0.0;
		cflux_zero = 0.0;
		ncont = 0.0;
		ncont_zero = 0.0;
		ncont_zero_scaled = 0.0;
		nflux = 0.0;
		nflux_zero = 0.0;
		pcont = 0.0;
		pcont_zero = 0.0;
		pcont_zero_scaled = 0.0;
		pflux = 0.0;
		pflux_zero = 0.0;
	}

	void init(Gridcell& gridcell);
	void check(Gridcell& gridcell);
	// indiv and patch-level functions are for use with true crop stands only
	void init_indiv(Individual& indiv);
	bool check_indiv(Individual& indiv, bool check_harvest = false);
	bool check_indiv_C(Individual& indiv, bool check_harvest = false);
	bool check_indiv_N(Individual& indiv, bool check_harvest = false);
	bool check_indiv_P(Individual& indiv, bool check_harvest = false);
	
	void init_patch(Patch& patch);
	bool check_patch(Patch& patch, bool check_harvest = false);
	bool check_patch_C(Patch& patch, bool check_harvest = false);
	bool check_patch_N(Patch& patch, bool check_harvest = false);
	bool check_patch_P(Patch& patch, bool check_harvest = false);

	void check_year(Gridcell& gridcell); // calls both check_year_C, check_year_N and check_year_P
	void check_year_C(Gridcell& gridcell);
	void check_year_N(Gridcell& gridcell);
	void check_year_P(Gridcell& gridcell);
	void check_period(Gridcell& gridcell);

	void serialize(ArchiveStream& arch);
};

/// This struct contains the result of a photosynthesis calculation.
/** \see photosynthesis */
struct PhotosynthesisResult : public Serializable {

	/// Constructs an empty result
	PhotosynthesisResult() {
		clear();
	}

	/// Clears all members
	/** This is returned by the photosynthesis function when no photosynthesis
	 *  takes place.
	 */
	void clear() {
		agd_g       = 0;
		adtmm       = 0;
		rd_g        = 0;
		vm          = 0;
		je          = 0;
		nactive_opt = 0.0;
		vmaxnlim    = 1.0;
		vmaxplim	= 1.0;
	}

	/// RuBisCO capacity (gC/m2/day)
	double vm;

	/// gross daily photosynthesis (gC/m2/day)
	double agd_g;

	/// leaf-level net daytime photosynthesis
	/** expressed in CO2 diffusion units (mm/m2/day) */
	double adtmm;

	/// leaf respiration (gC/m2/day)
	double rd_g;

	/// PAR-limited photosynthesis rate (gC/m2/h)
	double je;

	/// optimal leaf nitrogen associated with photosynthesis (kgN/m2)
	double nactive_opt;

	/// optimal leaf phosphorus associated with photosynthesis (kgP/m2)
	double pactive_opt;

	/// nitrogen limitation on vm
	double vmaxnlim;

	/// phosphorus limitation on vm
	double vmaxplim;

	/// net C-assimilation (gross photosynthesis minus leaf respiration) (kgC/m2/day)
    double net_assimilation() const {
		return (agd_g - rd_g) * KG_PER_G;
    }

	void serialize(ArchiveStream& arch);
};

/// Class containing serializable variables for Weathergenerator GWGen
/** 
Variables for the build-in random-generator and to keep track of whether past days 
were rain days
*/
class WeatherGenState : public Serializable {

public:
	/// Random state variable q
	int q[10];
	/// Random state variable carry
	int carry;
	/// Random state variable xcng
	int xcng;
	/// Random state variable xs
	unsigned int xs; 
	/// Random state variable indx
	int indx;
	/// Random state variable have
	bool have;
	/// Random state gamma	
	double gamma_vals[2];
	/// Indicator for whether the recent two days were rein-days
	bool pday[2];
	/// Random state's residuals
	double resid[4];

	WeatherGenState() {
		for (int i = 0; i<10; i++) q[i] = 0;
		carry = 0;
		xcng = 0;
		xs = 0;
		indx = 0;
		have = false;
		for (int i = 0; i<2; i++) gamma_vals[i] = 0.0;
		for (int i = 0; i<2; i++) pday[i] = false;
		for (int i = 0; i<4; i++) resid[i] = 0.0;
	}

	void serialize(ArchiveStream& arch);
};

/// This struct contains the environmental input to a photosynthesis calculation.
/** \see photosynthesis */
struct PhotosynthesisEnvironment {

	/// Constructs an empty result
	PhotosynthesisEnvironment() {
		clear();
	}

	/// Clears all members
	/** Nonsense values to cause a crash if used
	 */
	void clear() {
		co2 = 0;
		temp = 0;
		par = 0;
		fpar = 0;
		daylength = 0;
	}

private:
	/// atmospheric ambient CO2 concentration (ppmv)
	double co2;

	/// mean air temperature today (deg C)
	double temp;

	/// total daily photosynthetically-active radiation today (J / m2 / day) (ALPHAA not yet accounted for)
	double par;

	/// fraction of PAR absorbed by foliage
	double fpar;

	/// day length, must equal 24 in diurnal mode(h)
	double daylength;

public:
	/// set 
    void set(double co2_env, double temp_env, double par_env, double fpar_env, double daylength_env) {
		co2 = co2_env;
		temp = temp_env;
		par = par_env;
		fpar = fpar_env;
		daylength = daylength_env;
    }

	double get_apar() const {
		return par * fpar;
	}

	double get_par() const {
		return par;
	}

	double get_fpar() const {
		return fpar;
	}

	double get_temp() const {
		return temp;
	}

	double get_co2() const {
		return co2;
	}

	double get_daylength() const {
		return daylength;
	}

};

/// This struct contains the stresses used in a photosynthesis calculation.
/** \see photosynthesis */
struct PhotosynthesisStresses {

	/// Constructs an empty result
	PhotosynthesisStresses() {
		no_stress();
	}

	/// All members set to no stress values
	/** Default values indicating no stress
	 */
	void no_stress() {

		ifnlimvmax = false;
		ifplimvmax = false;
		moss_ps_limit = 1.0;
		graminoid_ps_limit = 1.0;
		inund_stress = 1.0;
	}

private:
	/// whether nitrogen should limit Vmax
	bool ifnlimvmax;

	/// whether phosphorus should limit Vmax
	bool ifplimvmax;

	///  limit to moss photosynthesis. [0,1], where 1 means no limit
	double moss_ps_limit;

	/// limit to graminoid photosynthesis. [0,1], where 1 means no limit
	double graminoid_ps_limit;

	/// limit to photosynthesis due to inundation, where 1 means no limit
	double inund_stress;

public:
	/// Set the stresses
    void set(bool thisnlimvmax, bool thisplimvmax, double this_moss_ps_limit, double this_graminoid_ps_limit, double this_inund_stress) {

		ifnlimvmax = thisnlimvmax;
		ifplimvmax = thisplimvmax;
		moss_ps_limit = this_moss_ps_limit;
		graminoid_ps_limit = this_graminoid_ps_limit;
		inund_stress = this_inund_stress;
    }

	bool get_ifnlimvmax() const {
		return ifnlimvmax;	
	}

	bool get_ifplimvmax() const {
		return ifplimvmax;
	}

	double get_moss_ps_limit() const {
		return moss_ps_limit;	
	}

	double get_graminoid_ps_limit() const {
		return graminoid_ps_limit;	
	}

	double get_inund_stress() const {
		return inund_stress;	
	}
};


/// The Climate for a grid cell
/** Stores all static and variable data relating to climate parameters, as well as
 *  latitude, atmospheric CO2 concentration and daylength for a grid cell. Includes
 *  a reference to the parent Gridcell object (defined below). Initialised by a
 *  call to initdrivers.
 */
class Climate : public Serializable {

	// MEMBER VARIABLES

public:
	/// reference to parent Gridcell object
	Gridcell& gridcell;

	/// values for randomisation in Weathergenerator GWGEN
	WeatherGenState weathergenstate;

	/// mean air temperature today (deg C)
	double temp;

	/// total daily net downward shortwave solar radiation today (J/m2/day)
	double rad;

	/// total daily photosynthetically-active radiation today (J/m2/day)
	double par;

	/// precipitation today (mm)
	double prec;

	/// 10 m wind (km/h)
	double u10;

	/// rel. humidity (fract.)
	double relhum;

	/// min and max daily temperature (deg C)
	double tmin, tmax; 

	/// day length today (h)
	double daylength;

	/// atmospheric ambient CO2 concentration today (ppmv)
	double co2;

	/// latitude (degrees; +=north, -=south)
	double lat;

	/// Insolation today, see also instype
	double insol;

	/// daily phosphorus deposition per month (kgP/m2/d)
	double mpdep[12];

	/// Type of insolation
	/** This decides how to interpret the variable insol,
	 *  see also documentation for the insoltype enum.
	 */
	insoltype instype;

	/// equilibrium evapotranspiration today (mm/day)
	double eet;

	/// mean temperature for the last 31 days (deg C)
	double mtemp;

	/// mean of lowest mean monthly temperature for the last 20 years (deg C)
	double mtemp_min20;

	/// mean of highest mean monthly temperature for the last 20 years (deg C)
	double mtemp_max20;

	/// highest mean monthly temperature for the last 12 months (deg C)
	double mtemp_max;

	/// accumulated growing degree day sum on 5 degree base
	/** reset when temperatures fall below 5 deg C */
	double gdd5;

	/// total gdd5 (accumulated) for this year (reset 1 January)
	double agdd5;

	/// accumulated growing degree day sum on 0 degree base (Wolf et al. 2008)
	double gdd0;

	/// total gdd0 (accumulated) for this year (reset 1 January)
	double agdd0;

	/// total gdd0 (accumulated) over each of the last 20 years
	/// climate.agdd0_20.mean() gives the average total gdd0 (accumulated) over the last 20 years
	Historic<double, 20> agdd0_20;

	/// number of days with temperatures <5 deg C
	/** reset when temperatures fall below 5 deg C;
	 *  maximum value is number of days in the year */
	int chilldays;

	/// true if chill day count may be reset by temperature fall below 5 deg C
	bool ifsensechill;

	/** Respiration response to today's air temperature incorporating damping of Q10
	 *  due to temperature acclimation (Lloyd & Taylor 1994)
	 */
	double gtemp;

	/// daily temperatures for the last 31 days (deg C)
	Historic<double, 31> dtemp_31;

	/// daily precipitation for the last 31 days (deg C)
	Historic<double, 31> dprec_31;

	/// daily eet for the last 31 days (deg C)
	Historic<double, 31> deet_31;

	/// minimum monthly temperatures for the last 20 years (deg C)
	double mtemp_min_20[20];

	/// maximum monthly temperatures for the last 20 years (deg C)
	double mtemp_max_20[20];

	/// minimum monthly temperature for the last 12 months (deg C)
	double mtemp_min;

	/// mean of monthly temperatures for the last 12 months (deg C)
	double atemp_mean;


	// BLAZE
	/// average annual rainfall (mm/a)
	double rainfall_annual_avg;
	/// current sum of annual Rainfall (mm)
	double rainfall_cur;
	/// Accumulated last rainfall (mm)
	double last_rainfall;
	/// Days since last rainfall 
	double days_since_last_rainfall;
	/// Keetch-Byram-Drought-Index
	double kbdi;
	/// McArthur forest fire index (FFDI)
	double mcarthur_forest_fire_index;	
	/// To keep track of running months daily FFDI 
	double ffdi_monthly[30];	

	// Saved parameters used by function daylengthinsoleet

	double sinelat;
	double cosinelat;
	double qo[Date::MAX_YEAR_LENGTH];
	double u[Date::MAX_YEAR_LENGTH];
	double v[Date::MAX_YEAR_LENGTH];
	double hh[Date::MAX_YEAR_LENGTH];
	double sinehh[Date::MAX_YEAR_LENGTH];
	double daylength_save[Date::MAX_YEAR_LENGTH];
	/// indicates whether saved values exist for this day
	bool doneday[Date::MAX_YEAR_LENGTH];

	/// diurnal temperature range, used in daily/monthly BVOC (deg C)
	double dtr;

	// containers for sub-daily values of temperature, short-wave downward
	// radiation, par, rad and gtemp (equivalent to temp, insol, par, rad and gtemp)
	// NB: units of these variable are the same as their daily counterparts,
	// i.e. representing daily averages (e.g. pars [J/m2/day])

	/// Sub-daily temperature (deg C) (\see temp)
	std::vector<double> temps;

	/// Sub-daily insolation (\see insol)
	std::vector<double> insols;

	/// Sub-daily PAR (\see par)
	std::vector<double> pars;

	/// Sub-daily net downward shortwave solar radiation (\see rad)
	std::vector<double> rads;

	/// Sub-daily respiration response (\see gtemp)
	std::vector<double> gtemps;

	/// Variables used for crop sowing date or seasonality calculation

	/// daily precipitations for the last 10 days (mm)
	double dprec_10[10];
	/// daily 10 day-sums of precipitations for today and yesterday (mm)
	double sprec_2[2];
	/// max temperature during the last test period
	double maxtemp;
	/// summer day when we test last year's crossing of sowing temperature limits; NH:June 30(day 180), SH:Dec.31(day 364), set in getgridcell()
	int testday_temp;
	/// last day of dry month when we test last year's crossing of sowing precipitation limits; NH:Dec.31(day 364), SH:June 30(day 180), set in getgridcell()
	int testday_prec;
	/// date used for sowing if no frost or spring occured during the year between the testmonths; NH:14, SH:195, set in getgridcell()
	int coldestday;
	/// used to adapt equations to hemisphere, set in getgridcell()
	int adjustlat;
	/// accumulated monthly pet values for this year
	double mpet_year[12];
	/// past 20 years monthly temperature values
	double mtemp_20[20][12];
	/// past 20 years monthly precipitation values
	double mprec_20[20][12];
	/// past 20 years monthly PET values
	double mpet_20[20][12];
	/// past 20 years monthly precipitation to PET ratios
	double mprec_pet_20[20][12];
	/// past 20 years minimum of monthly precipitation to PET ratios
	double mprec_petmin_20[20];
	/// past 20 years maximum of monthly precipitation to PET ratios
	double mprec_petmax_20[20];
	/// 20-year running average monthly temperature values
	double mtemp20[12];
	/// 20-year running average monthly precipitation values
	double mprec20[12];
	/// 20-year running average monthly PET values
	double mpet20[12];
	/// 20-year running average monthly precipitation to PET ratios
	double mprec_pet20[12];
	/// 20-year running average of minimum monthly precipitation to PET ratios
	double mprec_petmin20;
	/// 20-year running average of maximum monthly precipitation to PET ratios
	double mprec_petmax20;

	Historic<double, 20> hmtemp_20[12];
	Historic<double, 20> hmprec_20[12];
	Historic<double, 20> hmeet_20[12];

	/// seasonality type (SEASONALITY_NO, SEASONALITY_PREC, SEASONALITY_PRECTEMP, SEASONALITY_TEMP, SEASONALITY_TEMPPREC)
	seasonality_type seasonality;
	seasonality_type seasonality_lastyear;

	/// precipitation seasonality type (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
	/** based on the extremes of the 20-year monthly means
	 */
	prec_seasonality_type prec_seasonality;
	prec_seasonality_type prec_seasonality_lastyear;

	/// precipitation range (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
	/** based on the average of the 20-year monthly extremes
	 */
	prec_seasonality_type prec_range;
	prec_seasonality_type prec_range_lastyear;

	/// temperature seasonality (COLD, COLD_WARM, COLD_HOT, WARM, WARM_HOT, HOT)
	temp_seasonality_type temp_seasonality;
	temp_seasonality_type temp_seasonality_lastyear;

	/// whether several months with precipitation maxima exists (remains to be implemented)
	bool biseasonal;

	/// variation coefficient of 20-year mean monthly temperatures
	double var_prec;
	/// variation coefficient of 20-year mean monthly precipitation to PET ratios
	double var_temp;

	/// annual precipitation sum
	double aprec;
	/// annual average precipitation (last year) (mm)
	double aprec_lastyear;	

public:
	/// constructor function: initialises gridcell member
	Climate(Gridcell& gc):gridcell(gc) {

		var_prec = 0.0;
		var_temp = 0.0;

		aprec = 0.0;
		aprec_lastyear = 0.0;

		for(int m=0;m<12;m++) {

			mtemp20[m] = 0.0;
			mprec20[m] = 0.0;
			mpet20[m] = 0.0;
			mpet_year[m] = 0.0;
			mprec_pet20[m] = 0.0;

			for(int y=0;y<20;y++) {
				mtemp_20[y][m] = 0.0;
				mprec_20[y][m] = 0.0;
				mpet_20[y][m] = 0.0;
				mprec_pet_20[y][m] = 0.0;
			}
		}

		for(int y=0;y<20;y++) {

			mprec_petmin_20[y] = 0.0;
			mprec_petmax_20[y] = 0.0;
		}

		mprec_petmin20=0.0;
		mprec_petmax20=0.0;

		seasonality=SEASONALITY_NO;
		seasonality_lastyear=SEASONALITY_NO;
		prec_seasonality=DRY;
		prec_seasonality_lastyear=DRY;
		prec_range=DRY;
		prec_range_lastyear=DRY;
		temp_seasonality=COLD;
		temp_seasonality_lastyear=COLD;
		biseasonal=false;

		eet=0.0;

		gdd0 = 0.0;
		agdd0 = 0.0;
	};

	/// Initialises certain member variables
	/** Should be called before Climate object is applied to a new grid cell */
	void initdrivers(double latitude) {

		std::fill_n(mtemp_min_20, 20, 0.0);
		std::fill_n(mtemp_max_20, 20, 0.0);

		mtemp_min20 = 0.0;
		mtemp_max20 = 0.0;
		mtemp = 0.0;
		maxtemp = 0.0;
		gdd5 = 0.0;
		chilldays = 0;
		ifsensechill = true;
		atemp_mean = 0.0;

		lat = latitude;
		std::fill_n(doneday, Date::MAX_YEAR_LENGTH, false);
		sinelat = sin(lat * DEGTORAD);
		cosinelat = cos(lat * DEGTORAD);

		// Set crop-specific members
		if (latitude >= 0) {
			testday_temp = 180;		//June 30(day 180)
			testday_prec = 364;		//Dec.31(day 364)
			coldestday = COLDEST_DAY_NHEMISPHERE;
			adjustlat = 0;
		}
		else {
			testday_temp = 364;		//Dec.31(day 364)
			testday_prec = 180;		//June 30(day 180)
			coldestday = COLDEST_DAY_SHEMISPHERE;
			adjustlat = 181;
		}

		// BLAZE related variables
		rainfall_annual_avg      = 0.0; // average annual rainfall [mm]
		days_since_last_rainfall = 0;   // #Days-since-last-rainfall >3mm
		last_rainfall            = 0.0; // rainfall of last day of previous year [mm]
		kbdi                     = 0.0; // Keetch-Byram-Drought-index []

		for (int x=0; x<30; x++) {
			ffdi_monthly[x] = 0.;
		}

	}

	void serialize(ArchiveStream& arch);
};


/// Stores accumulated monthly and annual fluxes.
/** This class handles the storage and accounting of fluxes for a single patch.
 *  Different fluxes can be stored in different ways, depending on what kind of
 *  flux it is and what kind of output we want. The details of whether fluxes
 *  are stored per PFT or just as a patch total, or per day, month or only a
 *  yearly sum, is hidden from the 'scientific' code, which merely reports the
 *  fluxes generated.
 */
class Fluxes : public Serializable {

public:

	/// Fluxes stored as totals for the whole patch
	enum PerPatchFluxType {
		/// Carbon flux to atmosphere from burnt vegetation and litter (kgC/m2)
		FIREC,
		/// Carbon flux to atmosphere from soil respiration (kgC/m2)
		SOILC,
		/// Flux from atmosphere to vegetation associated with establishment (kgC/m2)
		ESTC,
		/// Flux to atmosphere from consumed harvested products (kgC/m2)
		HARVESTC,
		/// Flux from atmosphere to vegetation associated with sowing (kgC/m2)
		SEEDC,
		/// Flux from atmosphere to vegetation associated with manure addition (kgC/m2)
		MANUREC,
		/// Flux to vegetation associated with manure addition (kgN/m2)
		MANUREN,
		/// Flux to vegetation associated with N addition (kgN/m2) 
		NFERT,
		/// Nitrogen flux to atmosphere from consumed harvested products (kgN/m2)
		HARVESTN,
		/// Nitrogen flux from atmosphere to vegetation associated with sowing (kgN/m2)
		SEEDN,
		/// Flux to vegetation associated with manure addition (kgP/m2)
		MANUREP,
		/// Flux to vegetation associated with P addition (kgP/m2) 
		PFERT,
		/// Phosphorus flux to atmosphere from consumed harvested products (kgP/m2)
		HARVESTP,
		/// Phosphorus flux from atmosphere to vegetation associated with sowing (kgP/m2)
		SEEDP,
		/// NH3 flux to atmosphere from fire
		NH3_FIRE,
		/// NOx flux to atmosphere from fire
		NOx_FIRE,
		/// N2O flux to atmosphere from fire
		N2O_FIRE,
		/// N2 flux to atmosphere from fire
		N2_FIRE,
		/// P flux to atmosphere from fire
		P_FIRE,
		//---- Soil N transformation -----
		/// NH3 flux from soil (ntransform)
		NH3_SOIL,
		/// NO flux from soil (ntransform)
		NO_SOIL,
		/// N2O flux in soil (ntransform)
		N2O_SOIL,
		/// N2 flux from soil (ntransform)
		N2_SOIL,
		/// P flux from soil
		P_SOIL,
		/// DOC flux from soil (ntransform)
		DOC_FLUX,
		/// Net nitrification (ntransform)
		NET_NITRIF,
		/// Net denitrification (ntransform)
		NET_DENITRIF,
		/// Gross nitrification (ntransform)
		GROSS_NITRIF,
		/// Gross denitrification (ntransform)
		GROSS_DENITRIF,
		/// Reproduction costs
		REPRC,
		/// Total (i.e. CH4C_DIFF + CH4C_PLAN + CH4C_EBUL) CH4 flux to atmosphere from peatland soils (gC/m2).
		CH4C,
		/// Diffused CH4 flux to atmosphere from peatland soils (gC/m2).
		CH4C_DIFF,
		/// Plant-mediated CH4 flux to atmosphere from peatland soils (gC/m2).
		CH4C_PLAN,
		/// CH4 flux to atmosphere from peatland soils due to ebullition (gC/m2).
		CH4C_EBUL,
		/// Number of types, must be last
		NPERPATCHFLUXTYPES
	};

	/// Fluxes stored per pft
	enum PerPFTFluxType {
		/// NPP (kgC/m2)
		NPP,
		/// GPP (kgC/m2)
		GPP,
		/// Autotrophic respiration (kgC/m2)
		RA,
		/// Isoprene (mgC/m2)
		ISO,
		/// Monoterpene (mgC/m2)
		MT_APIN,
		MT_BPIN,
		MT_LIMO,
		MT_MYRC,
		MT_SABI,
		MT_CAMP,
		MT_TRIC,
		MT_TBOC,
		MT_OTHR,
		/// Number of types, must be last
		NPERPFTFLUXTYPES
	};

	// emission ratios from fire (NH3, NOx, N2O, N2) Delmas et al. 1995
	// values in .cpp file

	static const double NH3_FIRERATIO;
	static const double NOx_FIRERATIO;
	static const double N2O_FIRERATIO;
	static const double N2_FIRERATIO;

	/// Reference to patch to which this Fluxes object belongs
	Patch& patch;

	// MEMBER FUNCTIONS

public:
	/// constructor: initialises members
	Fluxes(Patch& p);

	/// Sets all fluxes to zero (call at the beginning of each year)
	void reset();

	void serialize(ArchiveStream& arch);

	/// Report flux for a certain flux type
	void report_flux(PerPFTFluxType flux_type, int pft_id, double value);

	/// Report flux for a certain flux type
	void report_flux(PerPatchFluxType flux_type, double value);

	/// \returns daily flux for a given flux type (for all PFTs)
	double get_daily_flux(PerPFTFluxType flux_type, int day) const;

	/// \returns daily flux for a given flux type
	double get_daily_flux(PerPatchFluxType flux_type, int day) const;

	/// \returns flux for a given month and flux type (for all PFTs)
	double get_monthly_flux(PerPFTFluxType flux_type, int month) const;

	/// \returns flux for a given month and flux type
	double get_monthly_flux(PerPatchFluxType flux_type, int month) const;

	/// \returns annual flux for a given PFT and flux type
	double get_annual_flux(PerPFTFluxType flux_type, int pft_id) const;

	/// \returns annual flux for a given flux type (for all PFTs)
	double get_annual_flux(PerPFTFluxType flux_type) const;

	/// \returns annual flux for a given flux type
	double get_annual_flux(PerPatchFluxType flux_type) const;

private:

	/// Stores one flux value per PFT and flux type
	std::vector<std::vector<double> > annual_fluxes_per_pft;

	/// Stores one flux value per month and flux type
	/** For the fluxes only stored as totals for the whole patch */
	double monthly_fluxes_patch[12][NPERPATCHFLUXTYPES];

	/// Stores one flux value per month and flux type
	/** For the fluxes stored per pft for annual values */
	double monthly_fluxes_pft[12][NPERPFTFLUXTYPES];

	/// Stores one flux value per day and flux type
	double daily_fluxes_patch[365][NPERPATCHFLUXTYPES];

	/// Stores one flux value per day and flux type
	double daily_fluxes_pft[365][NPERPFTFLUXTYPES];
};

/// Storage class of crop management information for one rotation period for a stand type, read from the instruction file.
class ManagementType {

public:
	/// id code (should be zero based and sequential, 0...nst-1)
	int id;
	/// name of management type
	xtring name;

	/// type of planting system ("", "MONOCULTURE", "SELECTION", etc.)
	xtring planting_system;
	/// type of harvest system ("", "CLEARCUT", "CONTINUOUS")
	xtring harvest_system;
	/// name of crop pft 
	xtring pftname;
	/// identifier of pft selection
	xtring selection;
	/// Rotation period in years
	double nyears;
	/// hydrology (RAINFED,IRRIGATED) 
	hydrologytype hydrology;
	/// irrigation efficiency
//	double firr;
	/// forced sowing date, unless sdate_force read from file
	int sdate;
	/// forced harvest date, unless hdate_force read from file
	int hdate;
	/// Nitrogen fertilisation amount, unless Nfert_read read from file
	double nfert;
	/// Phosphorus fertilisation amount, unless Pfert_read read from file
	double pfert;
	/// Whether grass is grown in fallow
	bool fallow;
	/// Double cropping of one crop (e.g. rice)
	bool multicrop;

	ManagementType() {

		planting_system = "";
		harvest_system = "";
		pftname = "";
		selection = "";
		nyears = 1.0;
		hydrology = RAINFED;
//		firr = 0.0;
		sdate = -1;
		hdate = -1;
		nfert = -1.0;
		pfert = -1.0;
		fallow = false;
		multicrop = false;
	}

	// Copy constructor
	ManagementType(const ManagementType& from) {

		name = from.name;
		pftname = from.pftname;
		hydrology = from.hydrology;
		sdate = from.sdate;
		hdate = from.hdate;
		nfert = from.nfert;
		pfert = from.pfert;
		fallow = from.fallow;
	}

	bool is_managed() {

		// Add new management parameters here
		if(pftname != "" || planting_system != "" || selection != ""||  harvest_system != "" ||  hydrology == IRRIGATED || fallow || nfert > -1.0 || pfert > -1.0)
			return true;
		else
			return false;
	}

	/// Returns true if pft is in pftselection.
	int pftinselection(const char* name) {

		bool found = false;
		char *p = NULL, string_copy[200] = {0};

		strcpy(string_copy, selection);
		p = strtok(string_copy, "\t\n ");
		if(p) {
			if(!strcmp(name, p)) {
				found = true;
			}
		}

		do {
			p = strtok(NULL, "\t\n ");
			if(p) {
				if(!strcmp(name, p)) {
					found = true;
				}
			}
		}
		while(p && !found);

		return found;
	}
};

/// A list of management types
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  management types objects.
 *
 * Functionality is inherited from the ListArray_id template type in the GUTIL
 * Library. Sequential management type objects can be referenced as array elements by id:
 *
 *   ManagementTypelist mtlist;
 *   ...
 *   for (i=0; i<nst; i++) {
 *     ManagementType& thismt=stlist[i];
 *     // query or modify object thismt here
 *   }
 *
 * or by iteration through the linked list:
 *
 *   mtlist.firstobj();
 *   while (mtlist.isobj) {
 *    ManagementType& thismt=mtlist.getobj();
 *     // query or modify object thismt here
 *     mtlist.nextobj();
 *   }
 */
class ManagementTypelist : public ListArray_id<ManagementType> {

public:
	int getmtid(xtring mtname) {

		int id = -1;

		for(unsigned int i=0; i< this->nobj; i++) {

			ManagementType& mt = (*this)[i];
			if(mt.name == mtname) {
				id = mt.id;
				break;
			}
		}

		return id;
	}
};

/// The one and only linked list of ManagementType objects	
extern ManagementTypelist mtlist;

/// Storage class of crop rotation information for a stand type, read from the instruction file.
struct CropRotation {

	/// Number of crops in rotation
	int ncrops;
	/// First rotation year
	int firstrotyear;

	CropRotation() {
		ncrops = 0;
		firstrotyear = 0;
	}
};

/// Stand type class for storing both static parameters, read from the instruction file,
/*	and dynamic variables, updated in landcover_change()
 *  Active stand types are stored in the stlist analogous to the pftlist.
 */
class StandType {

public:
	/// id code (should be zero based and sequential, 0...nst-1)
	int id;
	/// name of stand type
	xtring name;

	/// specifies type of landcover
	/** \see landcovertype */
	landcovertype landcover;	// specifies type of landcover (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
	/// Rotation information, read from the instruction file
	CropRotation rotation;
	/// Management struct (static)
	ManagementType management;
	/// Management types in a rotation cycle
	xtring mtnames[NROTATIONPERIODS_MAX];
	/// First management year: sets time when common features for managed stands begin, e.g. relaxed establishment rules and absence of disturbance before harvest begins
	/** \this currently only applies for stands with wood havest */
	int firstmanageyear;

	/// intercrop (NOINTERCROP,NATURALGRASS)
	intercroptype intercrop;
	/// whether natural pft:s are allowed to grow in stand type
	xtring naturalveg; // "", "GRASSONLY", "ALL"
	// whether only pft:s defined in management are allowed (plus intercrop or naturalveg/grass)
	bool restrictpfts;
	/// whether planted pft:s or all active pft:s are allowed to established after planting in a forest stand ("", "RESTRICTED", "ALL")
	xtring reestab;

	StandType() {

		intercrop = NOINTERCROP;
		naturalveg = "";
		restrictpfts = false;
		reestab = "ALL";
		firstmanageyear = 100000;
	}

	ManagementType& get_management(int rot = 0) {

		if(rotation.ncrops > 1) {
			return mtlist[mtlist.getmtid(mtnames[rot])];
		}
		else {
			return management;
		}
	}

	/// Returns position of management in rotation list if present. Returns -1 if not.
	int mtinrotation(xtring name) {

		int mtno = -1;
		for(int i=0; i<rotation.ncrops; i++) {
			if(name == mtnames[i])
				mtno = i;
		}

		return mtno;
	}

	/// Returns position of crop in rotation list if present. Returns -1 if not.
	int pftinrotation(xtring name) {

		int cropno = -1;
		for(int i=0; i<rotation.ncrops; i++) {
			if(name == get_management(i).pftname)
				cropno = i;
		}

		return cropno;
	}
};

/// A list of stand types
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  stand types objects.
 *
 * Functionality is inherited from the ListArray_id template type in the GUTIL
 * Library. Sequential stand type objects can be referenced as array elements by id:
 *
 *   StandTypelist stlist;
 *   ...
 *   for (i=0; i<nst; i++) {
 *     StandType& thisst=stlist[i];
 *     // query or modify object thisst here
 *   }
 *
 * or by iteration through the linked list:
 *
 *   stlist.firstobj();
 *   while (stlist.isobj) {
 *     StandType& thisst=stlist.getobj();
 *     // query or modify object thisst here
 *     stlist.nextobj();
 *   }
 */
class StandTypelist : public ListArray_id<StandType> {

};

/// The one and only linked list of StandType objects
extern StandTypelist stlist;


/// Holds static functional parameters for a plant functional type (PFT).
/** There should be one Pft object for each potentially occurring PFT. The same Pft object
 *  may be referenced (via the pft member of the Individual object; see below) by different
 *  average individuals. Member functions are included for initialising SLA given leaf
 *  longevity, and for initialising sapling/regen characteristics (required for
 * population mode).
 */
class Pft {

	// MEMBER VARIABLES

public:
	/// id code (should be zero based and sequential, 0...npft-1)
	int id;
	/// name of PFT
	xtring name;
	/// life form (tree or grass)
	lifeformtype lifeform;
	/// leaf phenology (raingreen, summergreen, evergreen, rain+summergreen, cropgreen)
	phenologytype phenology;
	/// leaf physiognomy (needleleaf, broadleaf)
	leafphysiognomytype leafphysiognomy;
	/// growing degree sum on 5 degree base required for full leaf cover
	double phengdd5ramp;
	/// water stress threshold for leaf abscission (range 0-1; raingreen PFTs)
	double wscal_min;
	/// biochemical pathway for photosynthesis (C3 or C4)
	pathwaytype pathway;
	/// approximate low temperature limit for photosynthesis (deg C)
	double pstemp_min;
	/// approximate lower range of temperature optimum for photosynthesis (deg C)
	double pstemp_low;
	/// approximate upper range of temperature optimum for photosynthesis (deg C)
	double pstemp_high;
	/// maximum temperature limit for photosynthesis (deg C)
	double pstemp_max;
	/// non-water-stressed ratio of intercellular to ambient CO2 partial pressure
	double lambda_max;
	/// vegetation root profile in an array containing fraction of roots in each soil layer, [0=upper layer]
	double rootdist[NSOILLAYER];
	/// shape parameter for initialisation of root distribtion
	double root_beta;
	/// canopy conductance component not associated with photosynthesis (mm/s)
	double gmin;
	/// maximum evapotranspiration rate (mm/day)
	double emax;
	/// maintenance respiration coefficient (0-1)
	double respcoeff;

	/// minimum leaf C:N mass ratio allowed when nitrogen demand is determined
	double cton_leaf_min;
	/// maximum leaf C:N mass ratio	allowed when nitrogen demand is determined
	double cton_leaf_max;
	/// average leaf C:N mass ratio (between min and max)
	double cton_leaf_avr;
	/// average fine root C:N mass ratio (connected cton_leaf_avr)
	double cton_root_avr;
	/// maximum fine root C:N mass ratio (used when mass is negligible)
	double cton_root_max;
	/// average sapwood C:N mass ratio (connected cton_leaf_avr)
	double cton_sap_avr;
	/// maximum sapwood C:N mass ratio (used when mass is negligible)
	double cton_sap_max;
	/// reference fine root C:N mass ratio
	double cton_root;
	/// reference sapwood C:N mass ratio
	double cton_sap;
	/// Maximum nitrogen (NH4+ and NO3- seperatly) uptake per fine root [kgN kgC-1 day-1]
	double nuptoroot;
	/// coefficient to compensate for vertical distribution of fine root on nitrogen uptake
	double nupscoeff;
	/// fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen longterm storage scalar
	double fnstorage;

	/// minimum leaf C:P mass ratio allowed when phosphorus demand is determined
	double ctop_leaf_min;
	/// maximum leaf C:P mass ratio	allowed when phosphorus demand is determined
	double ctop_leaf_max;
	/// average leaf C:P mass ratio (between min and max)
	double ctop_leaf_avr;
	/// average fine root C:P mass ratio (connected cton_leaf_avr)
	double ctop_root_avr;
	/// maximum fine root C:P mass ratio (used when mass is negligible)
	double ctop_root_max;
	/// average sapwood C:P mass ratio (connected cton_leaf_avr)
	double ctop_sap_avr;
	/// maximum sapwood C:P mass ratio (used when mass is negligible)
	double ctop_sap_max;
	/// reference fine root C:P mass ratio
	double ctop_root;
	/// reference sapwood C:P mass ratio
	double ctop_sap;
	/// Maximum phosphorus uptake per fine root [kgP kgC-1 day-1]
	double puptoroot;
	/// coefficient to compensate for vertical distribution of fine root on phosphorus uptake
	double pupscoeff;
	/// fraction of sapwood (root for herbaceous pfts) that can be used as a phosphorus longterm storage scalar
	double fpstorage;

	/// Michaelis-Menten kinetic parameters
	/** Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000) */
	double km_volume;
	/** Half saturation concentration for P uptake [kgP l-1] (Silveira & Cardoso 2004) */
	double kmp_volume;

	/// fraction of NPP allocated to reproduction
	double reprfrac;
	/// annual leaf turnover as a proportion of leaf C biomass
	double turnover_leaf;
	/// annual fine root turnover as a proportion of fine root C biomass
	double turnover_root;
	/// annual sapwood turnover as a proportion of sapwood C biomass
	double turnover_sap;
	/// sapwood and heartwood density (kgC/m3)
	double wooddens;
	/// maximum tree crown area (m2)
	double crownarea_max;
	/// constant in allometry equations
	double k_allom1;
	/// constant in allometry equations
	double k_allom2;
	/// constant in allometry equations
	double k_allom3;
	/// constant in allometry equations
	double k_rp;
	/// tree leaf to sapwood area ratio
	double k_latosa;
	/// specific leaf area (m2/kgC)
	double sla;
	/// leaf longevity (years)
	double leaflong;
	/// leaf to root mass ratio under non-water-stressed conditions
	double ltor_max;
	/// litter moisture flammability threshold (fraction of AWC)
	double litterme;
	/// fire resistance (0-1)
	double fireresist;
	/// minimum forest-floor PAR level for growth (grasses) or establishment (trees)
	/** J/m2/day, individual and cohort modes */
	double parff_min;
	/** parameter capturing non-linearity in recruitment rate relative to
	 *  understorey growing conditions for trees (Fulton 1991) (individual and
	 *  cohort modes)
	 */
	double alphar;
	/// maximum sapling establishment rate (saplings/m2/year) (individual and cohort modes)
	double est_max;
	/** constant used in calculation of sapling establishment rate when spatial
	 *  mass effect enabled (individual and cohort modes)
	 */
	double kest_repr;
	/// constant affecting amount of background establishment
	/** \see ifbgestab */
	double kest_bg;
	/** constant used in calculation of sapling establishment rate when spatial
	 *  mass effect disabled (individual and cohort modes)
	 */
	double kest_pres;
	/// expected longevity under non-stressed conditions (individual and cohort modes)
	double longevity;
	/// threshold growth efficiency for imposition of growth suppression mortality
	/** kgC/m2 leaf/year, individual and cohort modes */
	double greff_min;

	// Bioclimatic limits (all temperatures deg C)

	/// minimum 20-year coldest month mean temperature for survival
	double tcmin_surv;
	/// maximum 20-year coldest month mean temperature for establishment
	double tcmax_est;
	/// minimum degree day sum on 5 deg C base for establishment
	double gdd5min_est;
	/// minimum 20-year coldest month mean temperature for establishment
	double tcmin_est;
	/// minimum warmest month mean temperature for establishment
	double twmin_est;
	/// continentality parameter for boreal summergreen trees
	double twminusc;
	/// constant in equation for budburst chilling time requirement (Sykes et al 1996)
	double k_chilla;
	/// coefficient in equation for budburst chilling time requirement
	double k_chillb;
	/// exponent in equation for budburst chilling time requirement
	double k_chillk;
	/// array containing values for GDD0(c) given c=number of chill days
	/** Sykes et al 1996, Eqn 1
	 *  gdd0 has one element for each possible value for number of chill days
	 */
	double gdd0[Date::MAX_YEAR_LENGTH + 1];
	/// interception coefficient (unitless)
	double intc;

	/// the amount of N that is applied (kg N m-2)
	double N_appfert;
	/// the amount of P that is applied (kg P m-2)
	double P_appfert;
	/// 0 - 1 how much of the fertiliser is applied the first date, default 1.
	double fertrate[2];
	/// dates relative to sowing date
	int fertdates[2];
	double fert_stages[2];
	bool fertilised[2];

	double T_vn_min;
	double T_vn_opt;
	double T_vn_max;

	double T_veg_min;
	double T_veg_opt;
	double T_veg_max;

	double T_rep_min;
	double T_rep_opt;
	double T_rep_max;

	double photo[3];

	double dev_rate_veg;
	double dev_rate_rep;

	double a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3;
	double cton_stem_avr;
	double cton_stem_max;
	double ctop_stem_avr;
	double ctop_stem_max;

	/// Drought tolerance level (0 = very -> 1 = not at all) (unitless)
	/** Used to implement drought-limited establishment */
	double drought_tolerance;

	// bvoc

	/// aerodynamic conductance (m s-1)
	double ga;
	/// isoprene emission capacity (ug C g-1 h-1)
	double eps_iso;
	/// whether (1) or not (1) isoprene emissions show a seasonality
	bool seas_iso;
	/// monoterpene emission capacity (ug C g-1 h-1) per monoterpene species
	double eps_mon[NMTCOMPOUNDS];
	/// fraction of monoterpene production that goes into storage pool (-) per monoterpene species
	double storfrac_mon[NMTCOMPOUNDS];

	/// Bioclimatic limits parameters from Wolf et al. 2008

	/// snow max [mm]
	double max_snow;
	/// snow min [mm]
	double min_snow;
	/// GDD0 min
	double gdd0_min;
	/// GDD0 max
	double gdd0_max;
	
	/// New parameters from parameters from Wania et al. (2009a, 2009b, 2010)

	/// Days per month for which inundation is tolerated
	int inund_duration;
	/// Inundation stress is felt when the water table (mm) is above wtp_max
	double wtp_max;
	/// Whether this PFT has aerenchyma through which O2 and CH4 can be transported (Wania et al. 2010 - Sec 2.6)
	bool has_aerenchyma; 

	/// Sapling/regeneration characteristics (used only in population mode)
	/** For trees, on sapling individual basis (kgC); for grasses, on stand area basis,
	 *  kgC/m2 */
	struct {
		/// leaf C biomass
		double cmass_leaf;
		/// fine root C biomass
		double cmass_root;
		/// sapwood C biomass
		double cmass_sap;
		/// heartwood C biomass
		double cmass_heart;
	} regen;

	/// specifies type of landcover pft is allowed to grow in (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
	landcovertype landcover;
	/// pft selection
	xtring selection;
	/// fraction of residue outtake at harvest
	double res_outtake;
	/// harvest efficiency
	double harv_eff;
	/// harvest efficiency of intercrop grass
	double harv_eff_ic;
	/// fraction of harvested products that goes into patchpft.harvested_products_slow
	double harvest_slow_frac;
	/// yearly turnover fraction of patchpft.harvested_products_slow (goes to gridcell.acflux_harvest_slow)
	double turnover_harv_prod;
	/// whether pft may grow as cover crop
	bool isintercropgrass;
	/// whether autumn temperature dependent sowing date is calculated
	bool ifsdautumn;
	/// upper temperature limit for autumn sowing
	double tempautumn;
	/// lower temperature limit for spring sowing
	double tempspring;
	/// default length of growing period
	int lgp_def;
	/// upper minimum temperature limit for crop sowing
	double maxtemp_sowing;
	/// default sowing date in the northern hemisphere (julian day)
	int sdatenh;
	/// default sowing date in the southern hemisphere
	int sdatesh;
	/// whether sowing date adjusting equation is used
	bool sd_adjust;
	/// parameter 1 in sowing date adjusting equation
	double sd_adjust_par1;
	/// parameter 2 in sowing date adjusting equation
	double sd_adjust_par2;
	/// parameter 3 in sowing date adjusting equation
	double sd_adjust_par3;
	/// latest date for harvesting in the northern hemisphere
	int hlimitdatenh;
	/// latest date for harvesting in the southern hemisphere
	int hlimitdatesh;
	/// default base temperature (°C) for heat unit (hu) calculation
	double tb;
	/// temperature under which vernalisation is possible (°C)
	double trg;
	/// default number of vernalising days required
	int pvd;
    /// sensitivity to the photoperiod effect [0-1]
	double psens;
	/// basal photoperiod (h) (pb<ps for longer days plants)
	double pb;
	/// lag in days after sowing before vernalization starts
	int vern_lag;
	/// saturating photoperiod (h) (ps<pb for shorter days plants)
	double ps;
	/// default potential heat units required for crop maturity (degree-days)
	double phu;
	/// whether quadratic equation used for calculating potential heat units (Bondeau method)
	bool phu_calc_quad;
	/// whether linear equation used for calculating potential heat units (Bondeau method)
	bool phu_calc_lin;
	/// minimum potential heat units required for crop maturity (Bondeau method) (degree-days)
	double phu_min;
	/// maximum potential heat units required for crop maturity (Bondeau method) (degree-days)
	double phu_max;
	/// reduction factor of potential heat units in spring crops (Bondeau method) (degree-days)
	double phu_red_spring_sow;
	/// number of days of phu decrease in the linear phu equation (Bondeau method)
	double ndays_ramp_phu;
	/// intercept for the linear phu equation (Bondeau method)
	double phu_interc;
	/// fraction of growing season (phu) at which senescence starts [0-1]
	double fphusen;
	/// type of senescence curve (see Bondeau et al. 2007)
	bool shapesenescencenorm;
	/// fraction of maximal LAI still present at harvest [0-1]
	double flaimaxharvest;
	/// default maximum LAI (only used for intercrop grass in the case where no pasture is present in any stand)
	double laimax;
	/// whether harvestable organs are above ground
	bool aboveground_ho;
	/// optimum harvest index
	double hiopt;
	/// minimum harvest index
	double himin;
	/// initial fraction of growing season's npp allocated to roots
	double frootstart;
	/// final fraction of growing season's npp allocated to roots
	double frootend;
	/// autumn/spring sowing of pft:s with tempautumn = 1
	int forceautumnsowing;	//0 = NOFORCING,  1 = AUTUMNSOWING, 2 = SPRINGSOWING
	/// N limited version of pft
	bool nlim;
	/// P limited version of pft
	bool plim;

	double avg_cton(const double& min, const double& max) {
		return 2.0 / (1. / min + 1. / max);
	}

	double avg_ctop(const double& min, const double& max) {
		return 2.0 / (1. / min + 1. / max);
	}
	// MEMBER FUNCTIONS

public:

	/// Constructor (initialises array gdd0)
	Pft() {

		std::fill_n(gdd0, Date::MAX_YEAR_LENGTH + 1, -1.0); // value<0 signifies "unknown"; see function phenology()

		nlim = false;
		plim = false;
		root_beta = 0.0;

		drought_tolerance = 0.0; // Default, means that the PFT will never be limited by drought.
		res_outtake = 0.0;
		harv_eff = 0.0;
		harv_eff_ic = 0.0;
		turnover_harv_prod = 1.0;	// default 1 year turnover time

		isintercropgrass = false;
		ifsdautumn = false;
		maxtemp_sowing = 60;
		sdatenh = -1;
		sdatesh = -1;
		lgp_def = 190;
		hlimitdatenh = -1;
		hlimitdatesh = -1;
		tb = -999.9;
		trg = -999.9;
		pvd = -1;
		psens = -1.0;
		pb = -1.0;
		vern_lag=0;
		ps = -1.0;
		phu = -1.0;
		phu_red_spring_sow = 1.0;
		fphusen = -1.0;
		shapesenescencenorm = 0;
		flaimaxharvest = -1.0;
		laimax = 0.0;
		aboveground_ho = true;
		frootstart = 0.0;
		frootend = 0.0;
		forceautumnsowing = 0;

		fertrate[0] = 0.0;
		fertrate[1] = 1.0;
		fertdates[0] = 0;
		fertdates[1] = 30;

		fert_stages[0] = 0.5;
		fert_stages[1] = 0.9;
		fertilised[0] = false;
		fertilised[1] = false;

		N_appfert = 0.0;
		P_appfert = 0.0;

		T_vn_min=0.0;
		T_vn_opt=0.0;
		T_vn_max=0.0;
		T_veg_min=0.0;
		T_veg_opt=0.0;
		T_veg_max=0.0;
		T_rep_min=0.0;
		T_rep_opt=0.0;
		T_rep_max=0.0;

		a1=0.0;
		b1=0.0;
		c1=0.0;
		d1=0.0;
		a2=0.0;
		b2=0.0;
		c2=0.0;
		d2=0.0;
		a3=0.0;
		b3=0.0;
		c3=0.0;
		d3=0.0;

		for (int i=0; i<3; i++)
			photo[i]=0.0;
	}

	/// Calculates SLA given leaf longevity
	void initsla() {

		// SLA has to be supplied in the insfile for crops with N limitation
		if (!(phenology == CROPGREEN && (ifnlim || ifplim))) {

			// Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
			// m2/kgC)

			if (leafphysiognomy == BROADLEAF) {
				sla = 0.2 * pow(10.0, 2.41 - 0.38 * log10(12.0 * leaflong));
			}
			else if (leafphysiognomy == NEEDLELEAF) {
				sla = 0.2 * pow(10.0, 2.29 - 0.4 * log10(12.0 * leaflong));
			}
		}
	}

	/// Calculates minimum leaf C:N ratio given leaf longevity
	void init_cton_min() {
		// cton_leaf_min has to be supplied in the insfile for crops with N limitation
		if (!(phenology == CROPGREEN && (ifnlim || ifplim))) {
			// Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
			// kgN/kgC)

			if (leafphysiognomy == BROADLEAF)
				cton_leaf_min = 500.0 / pow(10.0, 1.75 - 0.33 * log10(12.0 * leaflong));
			else if (leafphysiognomy == NEEDLELEAF)
				cton_leaf_min = 500.0 / pow(10.0, 1.52 - 0.26 * log10(12.0 * leaflong));
		}
	}

	/// Calculates minimum leaf C:P ratio given leaf longevity (actually sla here, maybe change in the future, also for C:N)
	void init_ctop_min() {
		// ctop_leaf_min has to be supplied in the insfile for crops with P limitation
		if (!(phenology == CROPGREEN && (ifnlim || ifplim))) {
			// Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
			// kgN/kgC)

			if (leafphysiognomy == BROADLEAF)
				//ctop_leaf_min = Calculation from leaflong??;
				ctop_leaf_min = exp(8.63342 + log(sla) * -0.80936) / ((2.78 + 1.0) / 2.0); //Calculation from sla
			else if (leafphysiognomy == NEEDLELEAF)
				//ctop_leaf_min = Calculation from leaflong??;
				ctop_leaf_min = exp(8.63342 + log(sla) * -0.80936) / ((2.78 + 1.0) / 2.0); //Calculation from sla
		}
	}

	void init_cton_limits() {

		// Fraction between min and max C:N ratio White et al. 2000
		double frac_mintomax = (phenology == CROPGREEN && (ifnlim || ifplim)) ? 5.0 : 2.78;	// Use value also without nlim ?

		// Fraction between leaf and root C:N ratio
		double frac_leaftoroot = 1.16; // Friend et al. 1997

		// Fraction between leaf and sap wood C:N ratio
		double frac_leaftosap = 6.9;   // Friend et al. 1997

		// Max leaf C:N ratio
		cton_leaf_max = cton_leaf_min * frac_mintomax;

		// Average leaf C:N ratio
		cton_leaf_avr = avg_cton(cton_leaf_min, cton_leaf_max);

		// Tighter C:N ratio range for roots and sapwood: picked out thin air
		double frac_maxtomin = .9;

		// Maximum fine root C:N ratio
		cton_root_max = cton_leaf_max * frac_leaftoroot;

		double cton_root_min = cton_root_max * frac_maxtomin;

		// Average fine root C:N ratio
		cton_root_avr = avg_cton(cton_root_min, cton_root_max);

		// Maximum sap C:N ratio
		cton_sap_max  = cton_leaf_max * frac_leaftosap;

		double cton_sap_min = cton_sap_max * frac_maxtomin;

		// Average sap C:N ratio
		cton_sap_avr  = avg_cton(cton_sap_min, cton_sap_max);

		if (lifeform == GRASS || lifeform == MOSS) {
			respcoeff /= 2.0 * cton_root / (cton_root_avr + cton_root_min);
		} 
		else {
			respcoeff /= cton_root / (cton_root_avr + cton_root_min) +
			             cton_sap  / (cton_sap_avr  + cton_sap_min);
		}
		cton_stem_max = 1.0/(2.0*0.0034); //Maize params
		cton_stem_avr = 1.0/(2.0*0.0068);
	}

	void init_ctop_limits() {

		// Fraction between min and max C:N ratio White et al. 2000 SAME FOR P, MAKES SENSE?
		double frac_mintomax = (phenology == CROPGREEN && (ifnlim || ifplim)) ? 5.0 : 2.78;	// Use value also without nlim ?

																				// Fraction between leaf and root C:N ratio
		double frac_leaftoroot = 1.16; // Friend et al. 1997

									   // Fraction between leaf and sap wood C:N ratio
		double frac_leaftosap = 6.9;   // Friend et al. 1997

									   // Max leaf C:N ratio
		ctop_leaf_max = ctop_leaf_min * frac_mintomax;

		// Average leaf C:P ratio
		ctop_leaf_avr = avg_ctop(ctop_leaf_min, ctop_leaf_max);

		// Tighter C:P ratio range for roots and sapwood: picked out thin air
		double frac_maxtomin = .9;

		// Maximum fine root C:P ratio
		ctop_root_max = ctop_leaf_max * frac_leaftoroot;

		double ctop_root_min = ctop_root_max * frac_maxtomin;

		// Average fine root C:P ratio
		ctop_root_avr = avg_ctop(ctop_root_min, ctop_root_max);

		// Maximum sap C:P ratio
		ctop_sap_max = ctop_leaf_max * frac_leaftosap;

		double ctop_sap_min = ctop_sap_max * frac_maxtomin;

		// Average sap C:P ratio
		ctop_sap_avr = avg_ctop(ctop_sap_min, ctop_sap_max);

		if (lifeform == GRASS || lifeform == MOSS) {
			respcoeff /= 2.0 * ctop_root / (ctop_root_avr + ctop_root_min);
		}
		else {
			respcoeff /= ctop_root / (ctop_root_avr + ctop_root_min) +
				ctop_sap / (ctop_sap_avr + ctop_sap_min);
		}
		ctop_stem_max = 1.0 / (2.0*0.0034); //Maize params
		ctop_stem_avr = 1.0 / (2.0*0.0068);
	}

	/// Calculates coefficient to compensate for different vertical distribution of fine root on nitrogen uptake
	void init_nupscoeff() {

		// Fraction fine root in upper soil layer should have higher possibility for mineralized nitrogen uptake
		// Soil nitrogen profile is considered to have a exponential decline (Franzluebbers et al. 2009) giving
		// an approximate advantage of 2 of having more roots in the upper soil layer
		const double upper_adv = 2.0;

		// Simple solution until we get C and N in all soil layers.
		double rootdist_upper = 0.0;
		double rootdist_lower = 0.0;

		for (int sl = 0; sl < NSOILLAYER; sl++) {
			if (sl < NSOILLAYER_UPPER) 
				rootdist_upper += rootdist[sl];
			else
				rootdist_lower += rootdist[sl];
		}

		nupscoeff = rootdist_upper * upper_adv + rootdist_lower;

	}

	/// Calculates coefficient to compensate for different vertical distribution of fine root on phosphorus uptake
	void init_pupscoeff() {

		// Fraction fine root in upper soil layer should have higher possibility for mineralized nitrogen uptake
		// Soil nitrogen profile is considered to have a exponential decline (Franzluebbers et al. 2009) giving
		// an approximate advantage of 2 of having more roots in the upper soil layer
		const double upper_adv = 2.0;

		// Simple solution until we get C and N in all soil layers.
		double rootdist_upper = 0.0;
		double rootdist_lower = 0.0;

		for (int sl = 0; sl < NSOILLAYER; sl++) {
			if (sl < NSOILLAYER_UPPER)
				rootdist_upper += rootdist[sl];
			else
				rootdist_lower += rootdist[sl];
		}

		pupscoeff = rootdist_upper * upper_adv + rootdist_lower;

	}

	/// Initialises sapling/regen characteristics in population mode following LPJF formulation
	void initregen() {

		// see function allometry in growth module.

		// Note: primary PFT parameters, including SLA, must be set before this
		//       function is called

		const double REGENLAI_TREE = 1.5;
		const double REGENLAI_GRASS = 0.001;
		const double SAPLINGHW = 0.2;

		if (lifeform == TREE) {

			// Tree sapling characteristics

			regen.cmass_leaf = pow(REGENLAI_TREE * k_allom1 * pow(1.0 + SAPLINGHW, k_rp) *
				pow(4.0 * sla / PI / k_latosa, k_rp * 0.5) / sla, 2.0 / (2.0 - k_rp));

			regen.cmass_sap = wooddens * k_allom2 * pow((1.0 + SAPLINGHW) *
				sqrt(4.0 * regen.cmass_leaf * sla / PI / k_latosa), k_allom3) *
				regen.cmass_leaf * sla / k_latosa;

			regen.cmass_heart = SAPLINGHW * regen.cmass_sap;
		}
		else if (lifeform == GRASS || lifeform==MOSS) {

			// Grass regeneration characteristics

			regen.cmass_leaf = REGENLAI_GRASS / sla;
		}

			regen.cmass_root = 1.0 / ltor_max * regen.cmass_leaf;
	}
    
    // Inits root fractions in each soil layer through a shape parameter beta (see Jackson et al., 1996)
    void init_rootdist() {

		double depth = Dz_soil * CM_PER_MM;

		rootdist[0] = 1.0 - pow(root_beta, depth); // init first layer

        double tot = rootdist[0];
        for (int i=1; i<NSOILLAYER; i++){
			depth += Dz_soil * CM_PER_MM;
			rootdist[i] = 1.0 - pow(root_beta, depth) - (1.0 - pow(root_beta, depth - Dz_soil * CM_PER_MM));
            tot += rootdist[i];
        }
        // Calibrated the root_beta for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.
        // Sometimes the rootdist goes below our maximum soildepth. When that happens, put the residual fraction in lowest soil layer
        rootdist[NSOILLAYER-1] += 1.0 - tot;
    }
    
	bool ismoss() const {

		if (lifeform == MOSS)
			return true;
		else
			return false;
	}

	bool isgrass() const {

		if (lifeform == GRASS)
			return true;
		else
			return false;
	}

	bool istree() const {

		if (lifeform == TREE)
			return true;
		else
			return false;
	}

	bool iswetlandspecies() const {

		if (ismoss() || has_aerenchyma)
			return true;
		else
			return false;
	}
};


/// A list of PFTs
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  Pft objects. In general, there should be a single Pftlist object containing
 *  a single list of PFTs. Pft objects within the list are then referenced by the pft
 *  member of each Individual object.
 *
 * Functionality is inherited from the ListArray_id template type in the GUTIL
 * Library. Sequential Pft objects can be referenced as array elements by id:
 *
 *   Pftlist pftlist;
 *   ...
 *   for (i=0; i<npft; i++) {
 *     Pft& thispft=pftlist[i];
 *     // query or modify object thispft here
 *   }
 *
 * or by iteration through the linked list:
 *
 *   pftlist.firstobj();
 *   while (pftlist.isobj) {
 *     Pft& thispft=pftlist.getobj();
 *     // query or modify object thispft here
 *     pftlist.nextobj();
 *   }
 */
class Pftlist : public ListArray_id<Pft> {

public:
	int getpftid(xtring pftname) {

		int id = -1;

		for(unsigned int i=0; i< this->nobj; i++) {

			Pft& pft = (*this)[i];
			if(pft.name == pftname) {
				id = pft.id;
				break;
			}
		}

		return id;
	}
};

/// The one and only linked list of Pft objects
extern Pftlist pftlist;

/// Container for crop-specific data at the individual level
struct cropindiv_struct : public Serializable {

	//Plant carbon biomass variables are all on patch area basis (kgC/m2)

	/// year's harvestable organ C biomass (= ycmass_plant)
	double cmass_ho;
	/// above-ground pool C biomass (when calculating daily cmass_leaf from lai_crop) (= ycmass_agpool)
	double cmass_agpool;
	double cmass_stem;

	/// year's maximum value of leaf C biomass
	double cmass_leaf_max;
	/// grs_cmass_leaf value saved at day before senescence (for LAI-calculation in allometry)
	double cmass_leaf_sen;

	/// today's increase of whole plant C biomass
	double dcmass_plant;
	/// today's increase of leaf C biomass
	double dcmass_leaf;
	/// today's increase of root C biomass
	double dcmass_root;
	/// today's increase of harvestable organ C biomass
	double dcmass_ho;
	/// today's increase of above-ground pool C biomass
	double dcmass_agpool;
	double dcmass_stem;

	/// today's increase of leaf N biomass
	double dnmass_leaf;
	/// today's increase of root N biomass
	double dnmass_root;
	/// today's increase of harvestable organ N biomass
	double dnmass_ho;
	/// today's increase of above-ground pool N biomass
	double dnmass_agpool;

	/// today's increase of leaf P biomass
	double dpmass_leaf;
	/// today's increase of root P biomass
	double dpmass_root;
	/// today's increase of harvestable organ P biomass
	double dpmass_ho;
	/// today's increase of above-ground pool P biomass
	double dpmass_agpool;

	///CARBON
	/// daily updated whole plant C biomass, reset at harvest day
	double grs_cmass_plant;
	/// daily updated leaf C biomass, reset at harvest day
	double grs_cmass_leaf;
	/// daily updated root C biomass, reset at harvest day
	double grs_cmass_root;
	/// daily updated harvestable organ C biomass, reset at harvest day
	double grs_cmass_ho;
	/// daily updated above-ground pool C biomass, reset at harvest day
	double grs_cmass_agpool;
	/// daily updated dead leaf C biomass, reset at harvest day
	double grs_cmass_dead_leaf;
	/// daily updated stem pool C biomass, reset at harvest day
	double grs_cmass_stem;

	/// carbon content of harvestable organs saved on first day of land use change year
	double grs_cmass_leaf_luc;
	/// carbon content of harvestable organs saved on first day of land use change year
	double grs_cmass_root_luc;
	/// carbon content of harvestable organs saved on first day of land use change year
	double grs_cmass_ho_luc;
	/// carbon content of above-ground pool saved on first day of land use change year
	double grs_cmass_agpool_luc;
	/// carbon content of dead leaves saved on first day of land use change year
	double grs_cmass_dead_leaf_luc;
	/// carbon content of stem saved on first day of land use change year
	double grs_cmass_stem_luc;

	/// daily updated whole plant C biomass, reset at day 0
	double ycmass_plant;
	/// daily updated leaf C biomass, reset at day 0
	double ycmass_leaf;
	/// daily updated root C biomass, reset at day 0
	double ycmass_root;
	/// daily updated harvestable organ C biomass, reset at day 0
	double ycmass_ho;
	/// daily updated above-ground pool C biomass, reset at day 0
	double ycmass_agpool;
	/// daily updated dead leaf C biomass, reset at day 0
	double ycmass_dead_leaf;
	/// daily updated stem C biomass, reset at day 0
	double ycmass_stem;

	/// year's whole plant C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_plant;
	/// year's leaf C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_leaf;
	/// year's root C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_root;
	/// year's harvestable organ C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_ho;
	/// year's above-ground pool C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_agpool;
	/// year's stem C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_stem;

	///NITROGEN
	/// nitrogen content of harvestable organs
	double nmass_ho;
	/// nitrogen content of above-ground pool
	double nmass_agpool;
	/// nitrogen content of dead leaves
	double nmass_dead_leaf;

	/// nitrogen content of harvestable organs saved on first day of land use change year
	double nmass_ho_luc;
	/// nitrogen content of above-ground pool saved on first day of land use change year
	double nmass_agpool_luc;
	/// nitrogen content of dead leaves saved on first day of land use change year
	double nmass_dead_leaf_luc;

	/// daily updated leaf N biomass, reset at day 0
	double ynmass_leaf;
	/// daily updated root N biomass, reset at day 0
	double ynmass_root;
	/// daily updated harvestable organ N biomass, reset at day 0
	double ynmass_ho;
	/// daily updated above-ground pool N biomass, reset at day 0
	double ynmass_agpool;
	/// daily updated dead leaf N biomass, reset at day 0
	double ynmass_dead_leaf;

	/// year's leaf N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_leaf;
	/// year's root N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_root;
	/// year's harvestable organ N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_ho;
	/// year's above-ground pool N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_agpool;

	///PHOSPHORUS
	/// phosphorus content of harvestable organs
	double pmass_ho;
	/// nitrogen content of above-ground pool
	double pmass_agpool;
	/// nitrogen content of dead leaves
	double pmass_dead_leaf;

	/// nitrogen content of harvestable organs saved on first day of land use change year
	double pmass_ho_luc;
	/// nitrogen content of above-ground pool saved on first day of land use change year
	double pmass_agpool_luc;
	/// nitrogen content of dead leaves saved on first day of land use change year
	double pmass_dead_leaf_luc;

	/// daily updated leaf N biomass, reset at day 0
	double ypmass_leaf;
	/// daily updated root N biomass, reset at day 0
	double ypmass_root;
	/// daily updated harvestable organ N biomass, reset at day 0
	double ypmass_ho;
	/// daily updated above-ground pool N biomass, reset at day 0
	double ypmass_agpool;
	/// daily updated dead leaf N biomass, reset at day 0
	double ypmass_dead_leaf;

	/// year's leaf N biomass at time of harvest (cumulative if several harvest events)
	double harv_pmass_leaf;
	/// year's root N biomass at time of harvest (cumulative if several harvest events)
	double harv_pmass_root;
	/// year's harvestable organ N biomass at time of harvest (cumulative if several harvest events)
	double harv_pmass_ho;
	/// year's above-ground pool N biomass at time of harvest (cumulative if several harvest events)
	double harv_pmass_agpool;

	/// dry weight crop yield harvested this year (cumulative if several harvest events), based on harv_cmass_xx
	double harv_yield;

	/// harvestable organ C biomass at the last two harvest events this year
	double cmass_ho_harvest[2];
	/// harvestable organ N biomass at the last two harvest events this year
	double nmass_ho_harvest[2];
	/// harvestable organ P biomass at the last two harvest events this year
	double pmass_ho_harvest[2];
	/// dry weight crop yield at the last two harvest events this year
	double yield_harvest[2];

	/// dry weight crop yield grown this year (cumulative if several harvest events), based on ycmass_xx
	double yield;

	/// whether this pft is the main crop in the stand (pft.id==stand.pftid)
	bool isprimarycrop;
	/// whether this pft is allowed to compete with the main crop during the same growing period (for future use)
	bool isprimarycovegetation;
	/// whether this pft is grown during a second growing period, different from the primary (main) crop (for future use)
//	bool issecondarycrop;

	/// set to true if pft.isintercropgrass is true and the stand's main crop pft.intercrop is "naturalgrass"
	bool isintercropgrass;

	cropindiv_struct() {
		cmass_ho=0.0;
		cmass_agpool=0.0;
		cmass_stem = 0.0;
		cmass_leaf_max=0.0;
		cmass_leaf_sen=0.0;
		yield=0.0;
		yield_harvest[0]=0.0;
		yield_harvest[1]=0.0;
		dcmass_leaf=0.0;
		dcmass_root=0.0;
		dcmass_plant=0.0;
		dcmass_ho=0.0;
		dcmass_agpool=0.0;
		grs_cmass_leaf=0.0;
		grs_cmass_root=0.0;
		grs_cmass_plant=0.0;
		grs_cmass_ho=0.0;
		grs_cmass_agpool=0.0;
		grs_cmass_stem = 0.0;
		grs_cmass_dead_leaf = 0.0;
		grs_cmass_leaf_luc=0.0;
		grs_cmass_root_luc=0.0;
		grs_cmass_ho_luc=0.0;
		grs_cmass_agpool_luc=0.0;
		grs_cmass_dead_leaf_luc = 0.0;
		grs_cmass_stem_luc = 0.0;
		nmass_ho=0.0;
		nmass_agpool=0.0;
		nmass_dead_leaf = 0.0;
		pmass_ho = 0.0;
		pmass_agpool = 0.0;
		pmass_dead_leaf = 0.0;
		ycmass_leaf=0.0;
		ycmass_root=0.0;
		ycmass_plant=0.0;
		ycmass_ho=0.0;
		ycmass_agpool=0.0;
		ycmass_stem = 0.0;
		ycmass_dead_leaf = 0.0;
		harv_cmass_leaf=0.0;
		harv_cmass_root=0.0;
		harv_cmass_root=0.0;
		harv_cmass_ho=0.0;
		harv_yield=0.0;
		harv_cmass_agpool=0.0;
		harv_cmass_stem = 0.0;
		cmass_ho_harvest[0]=0.0;
		cmass_ho_harvest[1]=0.0;

		//Nitrogen
		dnmass_leaf=0.0;
		dnmass_root=0.0;
		dnmass_ho=0.0;
		dnmass_agpool=0.0;
		ynmass_leaf=0.0;
		ynmass_root=0.0;
		ynmass_ho=0.0;
		ynmass_agpool=0.0;
		ynmass_dead_leaf = 0.0;
		harv_nmass_leaf=0.0;
		harv_nmass_root=0.0;
		harv_nmass_root=0.0;
		harv_nmass_ho=0.0;
		harv_nmass_agpool=0.0;
		nmass_dead_leaf_luc = 0.0;
		nmass_ho_harvest[0]=0.0;
		nmass_ho_harvest[1]=0.0;

		//Phosphorus
		dpmass_leaf = 0.0;
		dpmass_root = 0.0;
		dpmass_ho = 0.0;
		dpmass_agpool = 0.0;
		ypmass_leaf = 0.0;
		ypmass_root = 0.0;
		ypmass_ho = 0.0;
		ypmass_agpool = 0.0;
		ypmass_dead_leaf = 0.0;
		harv_pmass_leaf = 0.0;
		harv_pmass_root = 0.0;
		harv_pmass_root = 0.0;
		harv_pmass_ho = 0.0;
		harv_pmass_agpool = 0.0;
		pmass_dead_leaf_luc = 0.0;
		pmass_ho_harvest[0] = 0.0;
		pmass_ho_harvest[1] = 0.0;

		isprimarycrop=false;
		isprimarycovegetation=false;
//		issecondarycrop=false;
		isintercropgrass=false;
	}

	void serialize(ArchiveStream& arch);
};


/// A vegetation individual.
/** In population mode this is the average individual of a PFT population;
 *  in cohort mode: the average individual of a cohort;
 *  in individual mode: an individual plant. Each grass PFT is represented as a single
 *  individual in all modes. Individual objects are collected within list arrays of
 *  class Vegetation (defined below), of which there is one for each patch, and include
 *  a reference to their 'parent' Vegetation object. Use the createobj member function
 *  of class Vegetation to add new individuals.
 */
class Individual : public Serializable {

public:
	/// reference to Pft object containing static parameters for this individual
	Pft& pft;
	/// reference to Vegetation object to which this Individual belongs
	Vegetation& vegetation;
	/// id code (0-based, sequential)
	int id;
	/// leaf C biomass on modelled area basis (kgC/m2)
	double cmass_leaf;
	/// fine root C biomass on modelled area basis (kgC/m2)
	double cmass_root;
	/// sapwood C biomass on modelled area basis (kgC/m2)
	double cmass_sap;
	/// heartwood C biomass on modelled area basis (kgC/m2)
	double cmass_heart;
	/// C "debt" (retrospective storage) (kgC/m2)
	double cmass_debt;
	/// Total C mass at land use change (kgC/m2)
	double cmass_tot_luc;
	/// leaf C mass after tunrnover
	double cmass_leaf_post_turnover;
	/// root C mass after turnover
	double cmass_root_post_turnover;
	/// Latest tunover day for this individual
	int last_turnover_day;

	/// leaf N biomass on modelled area basis (kgN/m2)
	double nmass_leaf;
	/// root N biomass on modelled area basis (kgN/m2)
	double nmass_root;
	/// sap N biomass on modelled area basis (kgN/m2)
	double nmass_sap;
	/// heart N biomass on modelled area basis (kgN/m2)
	double nmass_heart;

	/// leaf P biomass on modelled area basis (kgP/m2)
	double pmass_leaf;
	/// root P biomass on modelled area basis (kgN/m2)
	double pmass_root;
	/// sap P biomass on modelled area basis (kgN/m2)
	double pmass_sap;
	/// heart P biomass on modelled area basis (kgN/m2)
	double pmass_heart;

	/// leaf N biomass on modelled area basis saved on first day of land use change year
	double nmass_leaf_luc;
	/// root N biomass on modelled area basis on first day of land use change year
	double nmass_root_luc;
	/// sap N biomass on modelled area basis on first day of land use change year
	double nmass_sap_luc;
	/// heart N biomass on modelled area basis on first day of land use change year
	double nmass_heart_luc;
	/// total N biomass on modelled area basis on first day of land use change year
	double nmass_tot_luc;

	/// leaf P biomass on modelled area basis saved on first day of land use change year
	double pmass_leaf_luc;
	/// root P biomass on modelled area basis on first day of land use change year
	double pmass_root_luc;
	/// sap P biomass on modelled area basis on first day of land use change year
	double pmass_sap_luc;
	/// heart P biomass on modelled area basis on first day of land use change year
	double pmass_heart_luc;
	/// total P biomass on modelled area basis on first day of land use change year
	double pmass_tot_luc;

	/// foliar projective cover (FPC) under full leaf cover as fraction of modelled area
	double fpc;
	/// foliar projective cover (FPC) this day as fraction of modelled area
	double fpc_daily;
	/// fraction of PAR absorbed by foliage over projective area today, taking account of leaf phenological state
	double fpar;
	/// average density of individuals over patch (indiv/m2)
	double densindiv;
	/// vegetation phenological state (fraction of potential leaf cover)
	double phen;
	/// annual sum of daily fractional leaf cover
	/** Equivalent number of days with full leaf cover
	 *  (population mode only; reset on expected coldest day of year)
	 */
	double aphen;
	/// annual number of days with full leaf cover) (raingreen PFTs only; reset on 1 January)
	int aphen_raingreen;

	/// Photosynthesis values for this individual under non-water-stress conditions
	PhotosynthesisResult photosynthesis;

	/// sub-daily version of the above variable (NB: daily units)
	std::vector<PhotosynthesisResult> phots;

	/// accumulated NPP over modelled area (kgC/m2/year);
	/** annual NPP following call to growth module on last day of simulation year */
	double anpp;
	/// actual evapotranspiration over projected area (mm/day)
	double aet;
	/// annual actual evapotranspiration over projected area (mm/year)
	double aaet;
	/// leaf to root mass ratio
	double ltor;
	/// plant height (m)
	double height;
	/// plant crown area (m2)
	double crownarea;
	/// increment in fpc since last simulation year
	double deltafpc;
	/// bole height, i.e. height above ground of bottom of crown cylinder (m)
	/** (individual and cohort modes only) */
	double boleht;
	/// patch-level lai for this individual or cohort (function fpar)
	double lai;
	/// patch-level lai for cohort in current vertical layer (function fpar)
	double lai_layer;
	/// individual leaf area index (individual and cohort modes only)
	double lai_indiv;
	/// patch-level individual leaf area index (individual and cohort modes only)
	double lai_daily;
	/// daily individual leaf area index (individual and cohort modes only)
	double lai_indiv_daily;
	/// growth efficiency (NPP/leaf area) for each of the last five simulation years (kgC/m2/yr)
	Historic<double, NYEARGREFF> greff_5;
	/// increment of wood C for each of the last five simulation years (kgC/m2/yr)
	Historic<double, 10> cmass_wood_inc_5;
	/// individual/cohort age (years)
	double age;
	/// monthly LAI (including phenology component)
	double mlai[12];
	/// monthly maximum LAI (including phenology component)
	double mlai_max[12];
	/// FPAR assuming full leaf cover for all vegetation
	double fpar_leafon;
	/// LAI for current layer in canopy (cohort/individual mode; see function fpar)
	double lai_leafon_layer;
	/// non-water-stressed canopy conductance on FPC basis (mm/s)
	double gpterm;
	/// sub-daily version of the above variable (mm/s)
	std::vector<double> gpterms;
	/// interception associated with this individual today (patch basis)
	double intercep;

	/// accumulated mean fraction of potential leaf cover
	double phen_mean;

	/// whether individual subject to water stress
	bool wstress;

	/// leaf nitrogen that is photosyntetic active
	double nactive;
	/// Nitrogen extinction scalar
	/** Scalar to account for leaf nitrogen not following the optimal light
	  * extinction, but is shallower.
	  */
	double nextin;
	/// long-term storage of labile nitrogen
	double nstore_longterm;
	/// storage of labile nitrogen
	double nstore_labile;
	/// long-term storage of labile nitrogen saved on first day of land use change year
	double nstore_longterm_luc;
	/// storage of labile nitrogen saved on first day of land use change year
	double nstore_labile_luc;
	/// daily total nitrogen demand
	double ndemand;
	/// fraction of individual nitrogen demand available for uptake
	double fnuptake;
	/// annual nitrogen uptake
	double anuptake;
	/// maximum size of nitrogen storage
	double max_n_storage;
	/// scales annual npp to maximum nitrogen storage
	double scale_n_storage;
	/// annual nitrogen limitation on vmax
	double avmaxnlim;
	/// annual optimal leaf C:N ratio
	double cton_leaf_aopt;
	/// annual average leaf C:N ratio
	double cton_leaf_aavr;
	/// plant mobile nitrogen status
	double cton_status;
	/// total carbon in compartments before growth
	double cmass_veg;
	/// total nitrogen in compartments before growth
	double nmass_veg;
	/// whether individual optimal leaf nitrogen is above allowed limit
	bool n_opt_isabovelim;
	/// whether individual subject to nitrogen stress
	bool nstress;
	/// daily leaf nitrogen demand calculated from Vmax (kgN/m2)
	double leafndemand;
	/// daily root nitrogen demand based on leafndemand
	double rootndemand;
	/// daily sap wood nitrogen demand based on leafndemand
	double sapndemand;
	/// daily labile nitrogen demand based on npp
	double storendemand;
	/// daily harvestable organ nitrogen demand
	double hondemand;
	/// leaf fraction of total nitrogen demand
	double leaffndemand;
	/// root fraction of total nitrogen demand
	double rootfndemand;
	/// sap fraction of total nitrogen demand
	double sapfndemand;
	/// store fraction of total nitrogen demand
	double storefndemand;
	/// daily leaf nitrogen demand over possible uptake (storage demand)
	double leafndemand_store;
	/// daily root nitrogen demand over possible uptake (storage demand)
	double rootndemand_store;

	/// leaf phosphorus that is photosyntetic active
	double pactive;
	/// Phosphorus extinction scalar
	/** Scalar to account for leaf phosphorus not following the optimal light
	* extinction, but is shallower.
	*/
	double pextin;
	/// long-term storage of labile phosphorus
	double pstore_longterm;
	/// storage of labile phosphorus
	double pstore_labile;
	/// long-term storage of labile phosphorus saved on first day of land use change year
	double pstore_longterm_luc;
	/// storage of labile phosphorus saved on first day of land use change year
	double pstore_labile_luc;
	/// daily total phosphorus demand
	double pdemand;
	/// fraction of individual phosphorus demand available for uptake
	double fpuptake;
	/// annual phosphorus uptake
	double apuptake;
	/// maximum size of phosphorus storage
	double max_p_storage;
	/// scales annual npp to maximum phosphorus storage
	double scale_p_storage;
	/// annual phosphorus limitation on vmax
	double avmaxplim;
	/// annual optimal leaf C:P ratio
	double ctop_leaf_aopt;
	/// annual average leaf C:P ratio
	double ctop_leaf_aavr;
	/// plant mobile phosphorus status
	double ctop_status;
	/// total phosphorus in compartments before growth
	double pmass_veg;
	/// whether individual optimal leaf phosphorus is above allowed limit
	bool p_opt_isabovelim;
	/// whether individual subject to phosphorus stress
	bool pstress;
	/// daily leaf phosphorus demand calculated from Vmax (kgP/m2)
	double leafpdemand;
	/// daily root phosphorus demand based on leafpdemand
	double rootpdemand;
	/// daily sap wood phosphorus demand based on leafpdemand
	double sappdemand;
	/// daily labile phosphorus demand based on npp
	double storepdemand;
	/// daily harvestable organ phosphorus demand
	double hopdemand;
	/// leaf fraction of total phosphorus demand
	double leaffpdemand;
	/// root fraction of total phosphorus demand
	double rootfpdemand;
	/// sap fraction of total phosphorus demand
	double sapfpdemand;
	/// store fraction of total phosphorus demand
	double storefpdemand;
	/// daily leaf phosphorus demand over possible uptake (storage demand)
	double leafpdemand_store;
	/// daily root phosphorus demand over possible uptake (storage demand)
	double rootpdemand_store;

	/// The daily C lossed from leaves due to senescense, only crops.
	double daily_cmass_leafloss;
	/// The daily N lossed from leaves due to senescense, only crops.
	double daily_nmass_leafloss;
	/// The daily P lossed from leaves due to senescense, only crops.
	double daily_pmass_leafloss;
	/// The daily C lossed from roots due to senescense, only crops.
	double daily_cmass_rootloss;
	/// The daily N lossed from roots due to senescense, only crops.
	double daily_nmass_rootloss;
	/// The daily P lossed from roots due to senescense, only crops.
	double daily_pmass_rootloss;

	/// Number of days with non-negligible phenology this month
	int nday_leafon;
	// Whether this individual is truly alive.
	/** Set to false for first year after the Individual object is created, then true. */
	bool alive;
	/// NPP this day
	double dnpp;

	// bvoc

	/// isoprene production (mg C m-2 d-1)
	double iso;
	/// monoterpene production (mg C m-2 d-1) per monoteprene species
	double mon[NMTCOMPOUNDS];
	/// monoterpene storage pool (mg C m-2) per monoterpene species
	double monstor[NMTCOMPOUNDS];
	/// isoprene seasonality factor (-)
	double fvocseas;

	/// Pointer to struct with crop-specific data
	cropindiv_struct *cropindiv;

	// MEMBER FUNCTIONS

public:

	// Constructor function for objects of class Individual
	// Initialisation of certain member variables

	Individual(int i,Pft& p,Vegetation& v);
	~Individual();

	/// Access functions for cropindiv:
	cropindiv_struct* get_cropindiv() const;
	cropindiv_struct* set_cropindiv();

	void serialize(ArchiveStream& arch);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPFTFluxType flux_type, double value);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPatchFluxType flux_type, double value);

	/// Whether an individual is either a true crop or a cover crop grass
	inline bool istruecrop_or_intercropgrass() const {
		return (pft.landcover==CROPLAND && (pft.phenology==CROPGREEN || cropindiv->isintercropgrass));
	}

	/// Whether harvest and turnover is done on actual C and N on harvest or turnover day, which can occur any day of the year.
	bool has_daily_turnover() const;

	/// Whether resetting of grs_cmass and turnover (if has_daily_turnover() returns true) of continuous grass is to be done this day.
	/** This should occur at the very end of the growing period */
	bool is_turnover_day() const;

	/// Reduce current biomass due to mortality and/or fire
	/** The removed biomass is put into litter pools and/or goes to fire fluxes.
	 *
	 *  \param mortality      fraction of Individual's biomass killed due to
	 *                        mortality (including fire)
	 *  \param mortality_fire fraction of Individual's biomass killed due to
	 *                        fire only
	 */
	void reduce_biomass(double mortality, double mortality_fire);

	/// A version of the above reduce_biomass for the use with BLAZE
	void blaze_reduce_biomass(Patch& patch, double frac_survive);

	/// Total storage of nitrogen
	double nstore() const {
		return nstore_longterm + nstore_labile;
	}

	/// Total storage of phosphorus
	double pstore() const {
		return pstore_longterm + pstore_labile;
	}

	/// Total carbon wood biomass
	double cmass_wood() const {
		return cmass_sap + cmass_heart - cmass_debt;
	}

	/// Total nitrogen wood biomass
	double nmass_wood() const {
		return nmass_sap + nmass_heart;
	}

	/// Total phosphorus wood biomass
	double pmass_wood() const {
		return pmass_sap + pmass_heart;
	}

	/// Total carbon biomass
	double ccont(double scale_indiv = 1.0, bool luc = false) const;
	/// Total nitrogen biomass
	double ncont(double scale_indiv = 1.0, bool luc = false) const;
	/// Total phosphorus biomass
	double pcont(double scale_indiv = 1.0, bool luc = false) const;

	/// Whether grass growth is uninterrupted by crop growth.
	bool continous_grass() const;

	/// Checks whether any grs_cmass is negative, in which case it is zeroed and fluxes are corrected (only cropland).
	double check_C_mass();
	/// Checks whether any nmass is negative, in which case it is zeroed and fluxes are corrected (only cropland).
	double check_N_mass();
	/// Checks whether any nmass is negative, in which case it is zeroed and fluxes are corrected (only cropland).
	double check_P_mass();

	/// Save cmass-values on first day of the year of land cover change in expanding stands
	void save_cmass_luc();
	/// Save nmass-values on first day of the year of land cover change in expanding stands
	void save_nmass_luc();
	/// Save pmass-values on first day of the year of land cover change in expanding stands
	void save_pmass_luc();

	/// Current leaf C:N ratio
	/**
	 *  \param use_phen Set to false if indiv.phen shouldn't be considered
	 *                  when calculating C:N ratio
	 */
	double cton_leaf(bool use_phen = true) const;

	/// Current leaf C:P ratio
	/**
	*  \param use_phen Set to false if indiv.phen shouldn't be considered
	*                  when calculating C:P ratio
	*/
	double ctop_leaf(bool use_phen = true) const;

	/// Current fine root C:N ratio
	/**
	 *  \param use_phen Set to false if indiv.phen shouldn't be considered
	 *                  when calculating C:N ratio
	 */
	double cton_root(bool use_phen = true) const;

	/// Current fine root C:P ratio
	/**
	*  \param use_phen Set to false if indiv.phen shouldn't be considered
	*                  when calculating C:P ratio
	*/
	double ctop_root(bool use_phen = true) const;

	/// Current sap C:N ratio
	double cton_sap() const;

	/// Current sap C:P ratio
	double ctop_sap() const;

	/// Gets the individual's Patchpft
	Patchpft& patchpft() const;

	/// Transfers the individual's biomass (C and N) to litter and harvest pools/fluxes
	/**
	 *  \param harvest Set to true if some of the biomass should be harvested,
	 *                 harvest will be done according to the PFT's harvest efficiency
	 *                 and residue outtake.
	 */
	void kill(bool harvest = false);

	/// Annual mean wscal - water stress parameter (0-1 range; 1 = minimum stress)
	/** Value only valid at end of year, after call to canopy_exchange().
	 *
	 *  Currently, all Individuals belonging to a Patchpft share the same water stress.
	 */
	double wscal_mean() const;

	/// Gets the individual's daily cmass_leaf value
	double cmass_leaf_today() const;
	/// Gets the individual's daily cmass_root value
	double cmass_root_today() const;

	/// Gets the individual's daily LAI value (patch-level)
	/** Based on total leaf area for whatever the individual represents
	 *  (individual, cohort, population), over the whole patch.
	 */
	double lai_today() const;

	/// Gets the individual's daily LAI value (individual-level)
	/** Based on the leaf area for the average individual and
	 *  the average individual's crown area.
	 */
	double lai_indiv_today() const;

	/// Gets the Nitrogen limited LAI, Eq. 8 Olin 2015
	double lai_nitrogen_today() const;

	/// Gets the individual's daily fpc value
	double fpc_today() const;

	/// Gets the growingseason status for crop individual. Non-crop individuals always return true.
	bool growingseason() const;

	/// The N demand of the storage, only used for PNV.
	double ndemand_storage(double cton_leaf_opt);

	/// The P demand of the storage, only used for PNV.
	double pdemand_storage(double ctop_leaf_opt);
};


/// The vegetation in a patch - a list of individuals
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  Individual objects. A single Vegetation object is defined for each patch. A
 *  reference to the parent Patch object (defined below) is included as a member
 *  variable.
 *
 *  Functionality is inherited from the ListArray_idin1 template type in the GUTIL
 *  Library. Sequential Individual objects can be referenced as array elements by id,
 *  or by iteration through the linked list:
 *
 *    Vegetation vegetation
 *    ...
 *    vegetation.firstobj();
 *    while (vegetation.isobj) {
 *      Individual& thisindiv=vegetation.getobj();
 *      // query or modify object thisindiv here
 *      vegetation.nextobj();
 *    }
 */
class Vegetation : public ListArray_idin2<Individual,Pft,Vegetation>, public Serializable {

public:
	// MEMBER VARIABLES

	/// reference to parent Patch object
	Patch& patch;

	// MEMBER FUNCTIONS

	/// constructor (initialises member variable patch)
	Vegetation(Patch& p):patch(p) {};

	void serialize(ArchiveStream& arch);
};


/// Soiltype stores static parameters for soils and the snow pack.
/** One Soiltype object is defined for each Gridcell. State variables for soils
 *  are held by objects of class Soil, of which there is one for each patch
 *  (see below).
 */
class Soiltype {

	// MEMBER VARIABLES

public:

	/// fixed available water holding capacity of soil layers [0=upper layer] (mm)
	double awc[NSOILLAYER];
	/// fixed available water holding capacity of the standard Gerten soil layers [0=upper, 1 = lower] (mm)
	double gawc[2];

	/// coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996)
	double perc_base;
	/// coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996) for the evaporation soil layer
	double perc_base_evap;
	/// exponent in percolation calculation (=4 in Eqn 31, Haxeltine & Prentice 1996)
	double perc_exp;

	/// thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_0;
	/// thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_15;
	/// thermal diffusivity at 100% WHC (mm2/s)
	double thermdiff_100;

	/// wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	double wp[NSOILLAYER];
	/// saturation point. Cosby et al 1984
	double wsats[NSOILLAYER];
	
	/// organic soil fraction
	double org_frac_gridcell[NSOILLAYER];
	
	/// mineral soil fraction
	double min_frac_gridcell[NSOILLAYER];

	/// porosity of the soil
	double porosity_gridcell[NSOILLAYER];
	
	/// wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	// equivalents for the standard Gerten soil layers [0=upper, 1 = lower]
	double gwp[2];
	/// saturation point. Cosby et al 1984
	double gwsats[2];

	/// year at which to calculate equilibrium soil carbon
	int solvesom_end;
	/// year at which to begin documenting means for calculation of equilibrium soil carbon
	int solvesom_begin;

	/// water holding capacity plus wilting point for whole soil volume
	double wtot;

	// Sand, silt and clay fractions, should always add up to 1.
	/// fraction of soil that is sand
	double sand_frac;
	/// fraction of soil that is clay
	double clay_frac;
	/// fraction of soil that is silt
	double silt_frac;
	/// pH
	double pH;

	/// soilcode, 0 to 8
	int soilcode;
	/// volumetric fraction of organic material (m3 m-3) (Hillel, 1998)
	double organic_frac;
	// water held below wilting point, important for heat conductance
	// From the AGRMET Handbook, 2002
	double water_below_wp;
	/// porosity from AGRMET Handbook, 2002
	double porosity;
	// volumetric fraction of mineral material (m3 m-3) (Hillel, 1998)
	// = 1 - organic_frac - porosity
	double mineral_frac;
	// [cm], 3 depths at which monthly soil temperature is saved and output.
	// Defaults: 25, 75 and 150 cm (if array values == 0)
	double soiltempdepths[10];
	// Run on [mm/day]
	// Currently only used for wetlands, but could be used for irrigation too
	double runon;
	

	// PEAT PROPERTIES (needed if there is a peat stand in this gridcell)

	/// fixed available water holding capacity of soil layers [0=upper layer] (mm)
	double awc_peat[NSOILLAYER];
	/// fixed available water holding capacity of the standard Gerten soil layers [0=upper, 1 = lower] (mm)
	double gawc_peat[2];
	/// wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	double wp_peat[NSOILLAYER];
	/// saturation point. Cosby et al 1984
	double wsats_peat[NSOILLAYER];
	// equivalents for the standard Gerten soil layers [0=upper, 1 = lower]
	/// wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	double gwp_peat[2];
	/// saturation point. Cosby et al 1984
	double gwsats_peat[2];
	/// water holding capacity plus wilting point for whole soil volume
	double wtot_peat;
	/// fraction of soil that is sand
	double sand_frac_peat;
	/// fraction of soil that is clay
	double clay_frac_peat;
	/// fraction of soil that is silt
	double silt_frac_peat;
	///Soil Phosphorus parameters (moment fixed for Entisol)
	/// Empirical parameter for describing the equilibrium between labile P and sorbed P [kgP/m2]
	double kplab;
	/// maximum amount of sorbed P [kgP/m2]
	double spmax;
	/// annual Phosphorus weathering rate [kgP/m2/y]
	double pwtr;
	

	// MEMBER FUNCTIONS

public:

	/// Constructor: initialises certain member variables
	Soiltype() {

		solvesom_end = SOLVESOM_END;
		solvesom_begin = SOLVESOM_BEGIN;
		organic_frac = 0.02;
		pH = -1.0;
		
		// Assume no mineral content on peatlands
		sand_frac_peat = clay_frac_peat = silt_frac_peat = 0.0;

		runon = 0.0;
		for (int ii = 0; ii < 10; ii++)
			soiltempdepths[ii] = 0.0;
	}

	/// Override the default SOM years with 70-80% of the spin-up period length
	void updateSolveSOMvalues(const int& nyrspinup) {

		solvesom_end = static_cast<int>(0.8 * nyrspinup);
		solvesom_begin = static_cast<int>(0.7 * nyrspinup);

	}
};

/// CENTURY SOIL POOL
class Sompool : public Serializable {

public:

	/// Constructor
	Sompool() {

		// Initialise pool

		cmass = 0.0;
		nmass = 0.0;
		pmass = 0.0;
		ligcfrac = 0.0;
		delta_cmass = 0.0;
		delta_nmass = 0.0;
		delta_pmass = 0.0;
		fracremain = 0.0;
		litterme = 0.0;
		fireresist = 0.0;

		for (int m = 0; m < 12; m++) {
			mfracremain_mean[m] = 0.0;
		}
	}

	/// C mass in pool kgC/m2
	double cmass;
	/// Nitrogen mass in pool kgN/m2
	double nmass;
	/// Phosphorus mass in pool kgP/m2
	double pmass;
	/// (potential) decrease in C following decomposition today (kgC/m2)
	double cdec;
	/// (potential) decrease in nitrogen following decomposition today (kgN/m2)
	double ndec;
	/// (potential) decrease in phosphorus following decomposition today (kgP/m2)
	double pdec;
	/// daily change in carbon, nitrogen and phosphorus
	double delta_cmass, delta_nmass, delta_pmass;
	/// lignin fractions
	double ligcfrac;
	/// fraction of pool remaining after decomposition
	double fracremain;
	/// nitrogen to carbon ratio
	double ntoc;
	/// phosphorus to carbon ratio
	double ptoc;

	// Fire
	/// soil litter moisture flammability threshold (fraction of AWC)
	double litterme;
	/// soil litter fire resistance (0-1)
	double fireresist;

	// Fast SOM spinup variables

	/// monthly mean fraction of carbon pool remaining after decomposition
	double mfracremain_mean[12];

	void serialize(ArchiveStream& arch);
};

/// This struct contains litter for solving Century SOM pools.
/** \see equilsom() */
struct LitterSolveSOM : public Serializable {
	/// Constructs an empty result
	LitterSolveSOM() {
		clear();
	}

	/// Clears all members
	void clear() {
		for (int p = 0; p < NSOMPOOL; p++) {
			clitter[p] = 0.0;
			nlitter[p] = 0.0;
			plitter[p] = 0.0;
		}
	}

	/// Add litter
    void add_litter(double cvalue, double nvalue, double pvalue, int pool) {
		clitter[pool] += cvalue;
		nlitter[pool] += nvalue;
		nlitter[pool] += pvalue;
    }

	double get_clitter(int pool) {
		return clitter[pool];
	}
	double get_nlitter(int pool) {
		return nlitter[pool];
	}
	double get_plitter(int pool) {
		return plitter[pool];
	}

	void serialize(ArchiveStream& arch);

private:
	/// Carbon litter
	double clitter[NSOMPOOL];

	/// Nitrogen litter
	double nlitter[NSOMPOOL];

	/// Phosphorus litter
	double plitter[NSOMPOOL];
};

/// Soil stores state variables for soils and the snow pack.
/** Initialised by a call to initdrivers. One Soil object is defined for each patch.
 *  A reference to the parent Patch object (defined below) is included as a member
 *  variable. Soil static parameters are stored as objects of class Soiltype, of which
 *  there is one for each grid cell. A reference to the Soiltype object holding the
 *  static parameters for this soil is included as a member variable.
 *
 *  NB: The class Soil and its member functions and variables are declared in guess.h, 
 *      while its member functions are implemented in soil.cpp and in soilmethane.cpp.
 */
class Soil : public Serializable {

	// MEMBER VARIABLES
private:

	/// soil temperature today at 0.25 m depth (deg C)
	double temp25;
	/// water content of soil layers [0=upper layer] as fraction of available water holding capacity;
	double wcont[NSOILLAYER];
	/// water content of sublayer of upper soil layer for which evaporation from the bare soil surface is possible
	/** fraction of available water holding capacity */
	double wcont_evap;

public:
	/// reference to parent Patch object
	Patch& patch;
	/// reference to Soiltype object holding static parameters for this soil
	Soiltype& soiltype;
	/// the average wcont over the growing season, for each of the upper soil layers. Used in drought limited establishment. 
	double awcont_upper;
	/// daily water content in upper soil layer for each day of year
	double dwcontupper[Date::MAX_YEAR_LENGTH];
	/// mean water content in upper soil layer for last month
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double mwcontupper;
	/// stored snow as average over modelled area (mm rainfall equivalent)
	double snowpack;
	/// total runoff today (mm/day)
	double runoff;
	/// daily temperatures for the last month (deg C)
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double dtemp[31];
	/// mean soil temperature for the last month (deg C)
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double mtemp;
	/** respiration response to today's soil temperature at 0.25 m depth
	 *  incorporating damping of Q10 due to temperature acclimation (Lloyd & Taylor 1994)
	 */
	double gtemp;
	/// soil organic matter (SOM) pool with c. 1000 yr turnover (kgC/m2)
	double cpool_slow;
	/// soil organic matter (SOM) pool with c. 33 yr turnover (kgC/m2)
	double cpool_fast;

	// Running sums (converted to long term means) maintained by SOM dynamics module

	/// mean annual litter decomposition (kgC/m2/yr)
	double decomp_litter_mean;
	/// mean value of decay constant for fast SOM fraction
	double k_soilfast_mean;
	/// mean value of decay constant for slow SOM fraction
	double k_soilslow_mean;


	// Parameters used by function soiltemp and updated monthly

	double alag, exp_alag;

	/// water content of soil layers [0=upper layer] as fraction of available water holding capacity
	double mwcont[12][NSOILLAYER];
	/// daily water content in lower soil layer for each day of year
	double dwcontlower[Date::MAX_YEAR_LENGTH];
	/// mean water content in lower soil layer for last month
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double mwcontlower;

	/// rainfall and snowmelt today (mm)
	double rain_melt;
	/// upper limit for percolation (mm)
	double max_rain_melt;
	/// whether to percolate today
	bool percolate;

	//////////////////////////////////////////////////////////////////////////////////
	// Soil member variables

	/// sum of soil layers
	int ngroundl;
	/// density of the snowpack, daily
	double snowdens;
	

	// Temperature variables:
    
	/// temperature in each layer today (after call to cnstep) [deg C]
	double T[NLAYERS];
	/// Recorded T each day of the year, at each level [deg C]
	double T_soil[NLAYERS];
	/// Record the monthly average soil temp at SOILTEMPOUT layers [deg C]
	double T_soil_monthly[12][SOILTEMPOUT];
	/// soil temperature from previous time step
	double T_old[NLAYERS];
	/// soil temperature in each layer YESTERDAY
	double T_soil_yesterday[NLAYERS];
	/// soil temperature at 25 cm depth, as calculated using previous versions of the model [deg C]
	double temp_analyticsoln;
	

	// Padding - The values initialised in first call to calctemp method
    
	/// temperature in padding layers [deg C]
	double pad_temp[PAD_LAYERS];
	/// thickness of padding layers [mm]
	double pad_dz[PAD_LAYERS];
	

	// Ice and water variables for the soil layers, where Frac stands for Fraction
	
	/// (Frac)tion of ice in each layer: amount of ice / total volume of soil layer. Not associated with AWC. 
	double Frac_ice[NLAYERS];
	
	/// ice fraction in each layer YESTERDAY
	double Frac_ice_yesterday[NLAYERS];
	
	/// fraction of water in each layer: water / total volume of soil layer
	double Frac_water[NLAYERS];
	

	// Layer composition and properties:
	// Soil layer information:
    
	/// Thickness of soil layers [mm]
	double Dz[NLAYERS];
	/// porosity of each soil layer
	double por[NLAYERS];
	/// organic soil fraction
	double Frac_org[NLAYERS];
	/// peat fraction
	double Frac_peat[NLAYERS];
	/// mineral soil fraction
	double Frac_min[NLAYERS];
	/// Fraction of water held below permanent wilting point when no freexing
	double Fpwp_ref[NLAYERS];
	/// fraction of water held below permanent wilting point.
	double Frac_water_belowpwp[NLAYERS];
	/// diffusivity of the soil layers [mm2 / day]
	double Di[NLAYERS];
	/// heat capacities of the soil layers [J mm-3 K-1]
	double Ci[NLAYERS];
	/// thermal conductivities of the soil layers [J day-1 mm-1 K-1]
	double Ki[NLAYERS];

	bool snow_active; // active snow layer(s)?
	int snow_active_layers; // <= NLAYERS_SNOW
	/// liquid water in the snow pack [kg m-2] == [mm]
	double snow_water[NLAYERS_SNOW];
	/// ice in the snowpack [kg m-2]
	double snow_ice[NLAYERS_SNOW];

	// Log (l) of thermal conductivities. Used in the calculation of thermal conductivity and diffusivity.
	double lKorg, lKpeat, lKmin, lKwater, lKice, lKair;

	/// records the first time T is calculated
	bool firstTempCalc;


	// Daily storage:
    
	/// Depth in mm of the maximum depth of thaw this year
	double maxthawdepththisyear;
	/// daily thaw depth. The depth to the first soil layer with a temperature greater than 0 degrees C [mm]
	double thaw;
	/// monthly thawing depth full, where ALL the ice has melted [mm]
	double mthaw[12];
	/// daily thawing depth full, where ALL the ice has melted [mm]
	double dthaw[Date::MAX_YEAR_LENGTH];
	
	/// depth of the acrotelm [mm]
	double acro_depth;
	/// depth of the catotelm [mm]
	double cato_depth;
	/// acrotelm CO2 level [mimil L-1]
	double acro_co2;


	// Snow:
    
	/// days of continuous snow cover
	int snow_days;
	/// previous days of continuous snow cover
	int	snow_days_prev;
	/// daily snow depth [mm] 
	double dsnowdepth;
	/// Monthly snow depth (average) [mm]
	double msnowdepth[12];
	/// Previous December's snowdepth [mm] - used in establishment - from Wolf et al. (2008) 
	double dec_snowdepth;
  

	// Daily photosynthetic limits:
    
	/// Daily limit on moss photosynthetic activity due to dessication [0,1], but really [0.3,1]
	double dmoss_wtp_limit;
	/// Daily limit on graminoid photosynthetic activity as WTP drops.
	double dgraminoid_wtp_limit;
    /// acrotelm porosity MINUS Fgas (so, 0.98-0.08)
	double acro_por;
	/// acrotelm porosity MINUS Fgas (so, 0.92-0.08)
	double cato_por;
	

	// Peatland hydrology variables:
    
	/// daily water table position [mm]
	double wtp[Date::MAX_YEAR_LENGTH];
	/// monthly average water table position [mm]
	double mwtp[12];
	/// annual average water table position [mm]
	double awtp;
	/// water table depth from yesterday and updated today = -wtp [mm]
	double wtd;
	/// Water in acrotelm plus standing water (up to a max of 100mm) [mm]
	double Wtot;
	/// Standing water (up to a max of 100mm) [mm] - set to 0 in LPJG
	double stand_water;
	/// Volumetric water content in the NSUBLAYERS_ACRO of the acrotelm
	double sub_water[NSUBLAYERS_ACRO];
	// available water holding capacity of soil layers [0=upper layer] [mm], taking into
	// account the unavailability of frozen water. Default value: soiltype.awc[]
	double whc[NSOILLAYER];
	/// Max water (mm) that can be held in each layer
	double aw_max[NSOILLAYER];
	// Volumetric liquid water content. A fraction. Considers the entire (awc + Fpwp)
	// volumetric water content MINUS the ice fraction. Updated daily.
	double alwhc[NSOILLAYER];
	// Initial volumetric liquid water content. A fraction. Considers the entire
	// (awc + Fpwp) volumetric water content MINUS the ice fraction.
	double alwhc_init[NSOILLAYER];


	// Indices
    
	/// index for first active layer of soil
	int IDX;
	/// index for snow layers
	int SIDX;
	/// index for snow layer from previous day
	int SIDX_old;
	/// index for mixed layer
	int MIDX;
	
	// Hydrology variables:
    
	/// records the first time hydrology routine is called
	bool firstHydrologyCalc;
	
	// These are the number of sublayers in the standard 0.5/1.0m
	// hydrology laters. Set ONCE in hydrology routine
	int nsublayer1;
	int nsublayer2;
	int num_evaplayers;
	
	/// root fractions per layer
	double rootfrac[NLAYERS];
	/// air fraction in each layer
	double Frac_air[NLAYERS];
	
	/// daily carbon flux to atmosphere from soil respiration
	/// Temporary storage of heterotrophic respiration on PEATLAND until it is reduced by allocation of a certain fraction to CH4 production.
	double dcflux_soil;
	
	/// CH4 and CO2 stores in the soil layers - updated daily 
	double ch4_store;
	double co2_store;

	/// dissolved CO2 concentration in each layer [g CO2-C layer-1 d-1]
	double CO2_soil[NLAYERS];
	/// dissolved CO2 concentration in each layer yesterday [g CO2-C layer-1 d-1]
	double CO2_soil_yesterday[NLAYERS];
	/// daily CO2 production in each layer [g CO2-C layer-1 d-1]
	double CO2_soil_prod[NLAYERS];
	/// dissolved CH4 concentration in each layer [g CH4-C layer-1 d-1]
	double CH4[NLAYERS];
	/// dissolved CH4 concentration in each layer yesterday [g CH4-C layer-1 d-1]
	double CH4_yesterday[NLAYERS];
	/// daily CH4 production in each layer [g CH4-C layer-1 d-1]
	double CH4_prod[NLAYERS];
	/// daily CH4 oxidation in each layer [g CH4-C layer-1 d-1]
	double CH4_oxid[NLAYERS];
	/// CH4 which bubbles out [g CH4-C layer-1]
	double CH4_ebull_ind[NLAYERS];
	/// Volume of CH4 which bubbles out [m3]
	double CH4_ebull_vol[NLAYERS];
	/// gaseous CH4 concentration in each layer [g CH4-C layer-1 d-1]
	double CH4_gas[NLAYERS];
	/// dissolved CH4 concentration in each layer [g CH4-C layer-1 d-1]
	double CH4_diss[NLAYERS];
	/// gaseous CH4 concentration in each layer [g CH4-C layer-1] yesterday
	double CH4_gas_yesterday[NLAYERS];
	/// dissolved CH4 concentration in each layer [g CH4-C layer-1] yesterday
	double CH4_diss_yesterday[NLAYERS];
	/// gaseous CH4 volume in each layer [m3]
	double CH4_gas_vol[NLAYERS];
	// volumetric CH4 content [unitless]
	double CH4_vgc[NLAYERS];
	/// dissolved O2 concentration in each layer [mol O2 layer-1 d-1]
	double O2[NLAYERS];
	
	/// layer water volume [m3]
	double volume_liquid_water[NLAYERS];
	/// layer water + ice volume [m3]
	double total_volume_water[NLAYERS];
	/// tiller area
	double tiller_area[NLAYERS];
	

	// Gas diffusion variables:
    
	/// gas transport velocity of O2 [m d-1]
	double k_O2;
	/// gas transport velocity of CO2 [m d-1]
	double k_CO2;
	/// gas transport velocity of CH4 [m d-1]
	double k_CH4;
	/// Equilibrium concentration of O2 [mol L-1]
	double Ceq_O2;
	/// Equilibrium concentration of CO2 [mol L-1]
	double Ceq_CO2;
	/// Equilibrium concentration of CH4 [mol L-1]
	double Ceq_CH4;
	

//////////////////////////////////////////////////////////////////////////////////
// CENTURY SOM pools and other variables

	Sompool sompool[NSOMPOOL];

	/// daily percolation (mm)
	double dperc;
	/// fraction of decayed organic nitrogen leached each day;
	double orgleachfrac;

	/// soil NH4 mass in pool (kgN/m2)
	double NH4_mass;
	/// soil NO3 mass in pool (kgN/m2)
	double NO3_mass;
	/// soil NH4 mass input (kgN/m2)
	double NH4_input;
	/// soil NO3 mass input (kgN/m2)
	double NO3_input;
	/// annual sum of nitrogen mineralisation
	double anmin;
	/// annual sum of nitrogen immobilisation
	double animmob;
	/// annual leaching from available nitrogen pool
	double aminleach;
	/// annual leaching of organics from active nitrogen pool
	double aorgNleach;
	/// total annual nitrogen fixation
	double anfix;
	/// calculated annual mean nitrogen fixation
	double anfix_calc;
	/// annual leaching of organics nitrogen from carbon pool
	double aorgCleach;	

	/// soil PO4 mass in pool (kgP/m2)
	double pmass_labile;
	/// soil PO4 mass input (kgP/m2)
	double pmass_labile_input;
	/// soil PO4 mass delta for equilization with sorbed pool (kgP/m2)
	double pmass_labile_delta;
	/// annual sum of phosphorus mineralisation
	double apmin;
	/// annual sum of phosphorus immobilisation
	double apimmob;
	/// annual leaching from available phosphorus pool
	double aminpleach;
	/// annual leaching of organics from active phosphorus pool
	double aorgPleach;
	/// total annual phosphorus weathering
	double apwtr;
	/// total annual phosphorus deposition and fertilization
	double apdep;
	/// soil sorbed P pool  (kgP/m2)
	double pmass_sorbed;
	/// soil strongly sorbed P pool  (kgP/m2)
	double pmass_strongly_sorbed;
	/// soil occluded P pool  (kgP/m2)
	double pmass_occluded;

	// Variables for fast spinup of SOM pools

	/// monthly fraction of available mineral nitrogen taken up
	double fnuptake_mean[12];
	/// monthly fraction of organic carbon/nitrogen leached
	double morgleach_mean[12];
	/// monthly fraction of available mineral nitrogen leached
	double mminleach_mean[12];
	/// annual nitrogen fixation
	double anfix_mean;
	/// annual phosphorus weathering
	double apwtr_mean;
	/// annual phosphorus deposition and fertilization
	double apdep_mean;

	/// monthly fraction of available mineral phosphorus taken up
	double fpuptake_mean[12];
	/// monthly fraction of organic carbon/phosphorus leached
	double morgPleach_mean[12];
	/// monthly fraction of available mineral phosphorus leached
	double mminpleach_mean[12];

	// Solving Century SOM pools

	/// years at which to begin documenting for calculation of Century equilibrium
	int solvesomcent_beginyr;
	/// years at which to end documentation and start calculation of Century equilibrium
	int solvesomcent_endyr;

	/// Cumulative litter pools for one year.
	LitterSolveSOM litterSolveSOM;

	std::vector<LitterSolveSOM> solvesom;

	/// stored NH4 deposition in snowpack
	double snowpack_NH4_mass;
	/// stored NO3 deposition in snowpack
	double snowpack_NO3_mass;
	/// stored PO4 deposition in snowpack
	double snowpack_pmass_labile;

	/// pools of soil N species in transformation (nitrification & denitrifiacation)

	/// soil NH4 mass in pool (kgN/m2)
	// double NH4_mass;	// total, definde above
	double NH4_mass_w;	// wet proportion
	double NH4_mass_d;	// dry...

	/// soil NO3 mass in pool (kgN/m2)
	// double NO3_mass; // total, definde above
	double NO3_mass_w;
	double NO3_mass_d;
	/// soil NO2 mass in pool (kgN/m2)
	double NO2_mass;
	double NO2_mass_w;
	double NO2_mass_d;
	/// soil NO mass in pool (kgN/m2)
	double NO_mass;
	double NO_mass_w;
	double NO_mass_d;
	/// soil NO mass in pool (kgN/m2)
	double N2O_mass;
	double N2O_mass_w;
	double N2O_mass_d;
	/// soil N2 mass in pool (kgN/m2)
	double N2_mass;

	// soil pH
	double pH;	//TODO: pH - not used yet. Daily mean precip, based on annual average

	// soil labile carbon availability daily (kgC/m2/day)
	double labile_carbon;
	double labile_carbon_w;
	double labile_carbon_d;

	// MEMBER FUNCTIONS

public:
	/// constructor (initialises member variable patch)
	Soil(Patch& p,Soiltype& s):patch(p),soiltype(s) {
			init_states();
	}

	void init_states();

	/// return soil temperature at 25cm depth
	double get_soil_temp_25() const;

	/// Analytic soil temperature calculation for this patch. Updates daily
	void soil_temp_analytic(const Climate& climate, double depth);

	/// main hydrology routine
	void hydrology_lpjf(const Climate& climate, double fevap);	

	/// simpler hydrology routine from LPJ-GUESS v4.0
	void hydrology_lpjf_twolayer(const Climate& climate, double fevap);

	/// return wcont for a certain layer
	double get_layer_soil_water(int layer) const;

	/// return wcont_evap
	double get_layer_soil_water_evap() const;

	/// method to return the 'old' (i.e. Gerten) wcont for layers between layer1 and layer2
	double get_soil_water(int layer1, int layer2) const;

	/// method to return the 'old' (i.e. Gerten) wcont for upper layers
	double get_soil_water_upper() const;

	/// method to return the 'old' (i.e. Gerten) wcont for lower layers
	double get_soil_water_lower() const;

	/// copy wcont
	void copy_layer_soil_water_array(double wconttoreturn[NSOILLAYER]);

	/// add water to a certain soil layer and return the water that cannot be added due to limited capacity
	double add_layer_soil_water(int layer, double extrawater);

	/// set wcont for a certain layer
	void set_layer_soil_water(int layer, double newlayerwater);

	/// set wcont_evap
	void set_layer_soil_water_evap(double newevaplayerwater);

	/// update wcont_evap, whc[], Frac_water etc. based on wcont
	void update_soil_water();

	/// return true if there is more than 5% ice content in any ice in the top 50cm of soil (needed for irrigation)
	bool ice_in_top_layer();

	/// Peatland hydrology routine. Implements the peatland hydrology scheme of Wania et al. (2008)
	void hydrology_peat(const Climate& climate, double fevap);

	/// Main soil temperature calculation. See soil.cpp for definition
	bool soil_temp_multilayer(const double& temp);

	/// Calculate methane dynamics today
	bool methane(bool generatemethane);
	
	/// Initialise the root fractions in each layer of the wetland - see Wania et al. (2010)
	void init_peatland_root_fractions();

	/// Get soil clay fraction
	double get_clayfrac();

	/// Get soil sand fraction
	double get_sandfrac();

	/// Get soil silt fraction
	double get_siltfrac();

	/// Return CH4 content (g CH4-C / m2) 
	double get_ch4_content();

	/// Return CO2 content (gC / m2) 
	double get_co2_content();

	/// Water filled pore space of layer ([0-1])
	double wfps(int layer) const;

	/// Soil water freezing allowed?
	bool can_freeze() const;

	// Soil helper functions
	double nmass_avail(int pref = NO);
	void nmass_subtract(double nmass, int pref = NO);
	void nmass_inc(double nmass, int pref = NO);
	void nmass_multiplic_inc(double inc, int pref = NO);

	//// Soil helper functions
	//double pmass_avail(int pref = NO);
	//void pmass_subtract(double pmass, int pref = NO);
	//void pmass_inc(double pmass, int pref = NO);
	//void pmass_multiplic_inc(double inc, int pref = NO);

private:

	// Private member variables.
	static const bool ifallowphasechanges = true;
	static const bool snowdensityconstant = false;

	// Private helper methods. Contain daily functionality moved from calctemp. 
	// All definitions in soil.cpp and methane.cpp

	// Soil temperature and hydrology methods
	void update_from_yesterday();
	void update_snow_properties(const int& daynum, const double& dailyairtemp, double& Dsnow, double& Csnow, double& Ksnow);
	void snowpack_dynamics(const double &snowdepth, const int& soilsurfaceindex, int& snow_active_layers);
	void update_soil_diffusivities(const int& daynum, bool ansoln);
	void update_ice_fraction(const int& daynum, const int& MIDX);
	void update_layer_fractions(const int& daynum, const int& mixedl, const int& MIDX);

	// Hydrology methods
	void init_hydrology_variables();
	bool update_layer_water_content(int day);
	void update_acrotelm_co2(double atmo_co2);
	bool valid_layer_num(const int& ngr) const;

	// CH4 methods
	void calculate_carbon_store(int dy, bool today);
	bool calculate_gas_diffusivities(double CH4[NLAYERS], double CO2[NLAYERS], double O2[NLAYERS]);
	bool update_daily_gas_parameters();
	double diffuse_gas(double Cgas[NLAYERS], double D[NLAYERS], gastype thisgastype, double Ceq,
		double kgas, double Dz_m[NLAYERS], double &dailyDiff);
	bool plant_gas_transport(double Cgas[NLAYERS], double Ceq, double kgas, gastype thisgastyp,
		double& plantCH4TransportToday);
	double calculate_tiller_areas(double r_frac[NLAYERS], double t_area[NLAYERS]);
	bool calculate_gas_ebullition(double& ebull_today);

	void serialize(ArchiveStream& arch);
};

/// Container for crop-specific data at patchpft level
struct cropphen_struct : public Serializable {

	/// latest sowing date
	int sdate;
	/// sowing date of growing period ending in latest harvest this year
	int sdate_harv;
	/// sowing dates of growing periods ending in the two latest harvests this year
	int sdate_harvest[2];
	/// sowing dates of growing periods starting this year
	int sdate_thisyear[2];
	/// number of sowings this year
	int nsow;
	/// latest harvest date
	int hdate;
	/// two latest harvest dates this year
	int hdate_harvest[2];
	/// last date for harvest
	int hlimitdate;
	/// last day of heat unit sampling period, set in Crop_sowing_date_new()
	int hucountend;
	/// number of harvests this year
	int nharv;
	/// whether sdate_harvest[0] happened last year
	bool sownlastyear;
	/// latest senescence start date this year
	int sendate;
	/// latest beginning of intercropseason (2 weeks after the harvest date)
	int bicdate;
	/// latest end of intercropseason (2 weeks before the sowing date)
	int eicdate;
	/// number of growing days this growing period
	int growingdays;
	/// number of growing days this year (used for wscal_mean calculation)
	int growingdays_y;
	/// length of growingseason ending in last harvest
	int lgp;
	/// base temp for heat unit calculation (°C)
	double tb;
	/// number of vernalising days required
	int pvd;
	/// number of accumulated vernalizing days
	int vdsum;
	/// heat unit reduction factor due to vernalization [0-1]
	double vrf;
	/// heat unit reduction factor due to photoperiodism [0-1]
	double prf;
	/// potential heat units required for crop maturity (°Cd)
	double phu;
	/// potential heat units that would have been used without dynamic phu calculation
	double phu_old;
	/// heat unit sum aquired during last growing period (°Cd)
	double husum;
	/// heat unit sum aquired durin sampling period, starting with sdate
	double husum_sampled;
	/// this year's heat unit sum aquired from sdate to hucountend
	double husum_max;
	/// running mean of recent past's husum_max
	double husum_max_10;
	/// number of heat units sampling years
	int nyears_hu_sample;
	/// fraction of growing season [0-1] (husum/phu)
	double fphu;
	/// fraction of growing season at latest harvest
	double fphu_harv;
	/// whether in period of heat unit sampling
	bool hu_samplingperiod;
	/// number of heat unit sampling days
	int hu_samplingdays;
	/// harvest index today [0-1, >1 if below-ground ho], harvestable organ/above-ground C for above-ground harvestable organs, dependent on fphu, reduced by water stress
	double hi;
	/// fraction of harvest index today
	double fhi;
	/// phenology (fphu) contribution of fraction of harvest index today
	double fhi_phen;	//Phenology (fPHU) compoment of fhi
	/// water stress contribution of fraction of harvest index today
	double fhi_water;
	/// fraction of harvest index at latest harvest
	double fhi_harv;
	/// sum of crop patch demand (patch.wdemand) during crop growing period, reset on harvest day
	double demandsum_crop;
	/// sum of crop supply (patchpft.wsupply) during crop growing period, reset on harvest day
	double supplysum_crop;

	/// whether inside crop/intercrop grass growing period
	bool growingseason;
	/// whether yesterday was inside crop/intercrop grass growing period
	bool growingseason_ystd;
	/// whether inside crop senescence
	bool senescence;
	/// whether yesterday was inside crop senescence
	bool senescence_ystd;
	/// whether inside intercrop crass growing period (main crop pft variable)
	bool intercropseason;

	double vdsum_alloc;
	double vd;

	/// The fraction of the daily assimilates allocated to roots.
	double f_alloc_root;
	/// The fraction of the daily assimilates allocated to leaves.
	double f_alloc_leaf;
	/// The fraction of the daily assimilates allocated to harvestable organs, seeds.
	double f_alloc_horg;
	/// The fraction of the daily assimilates allocated to stem.
	double f_alloc_stem;
	/// Development stage from Wang & Engel 1998
	double dev_stage;
	// A variable holding the memory of whether this field was fertilised or not.
	bool fertilised[3];

	cropphen_struct() {
		sdate=-1;
		sdate_harv=-1;
		nsow=0;
		sownlastyear=false;
		sendate=-1;
		hdate=-1;
		hlimitdate=-1;
		hucountend=-1;
		nharv=0;
		tb=0.0;
		pvd=0;
		vdsum=0;
		vrf=1.0;
		phu=0.0;
		phu_old=0.0;
		husum_max=0.0;
		husum_sampled=0.0;
		husum_max_10=0.0;
		nyears_hu_sample = 0;
		prf=1.0;
		husum=0.0;
		fphu=0.0;
		fphu_harv=0.0;
		hu_samplingdays=0;
		hu_samplingperiod=false;

		hi=0.0;
		fhi=0.0;
		fhi_phen=0.0;
		fhi_water=1.0;
		fhi_harv=0.0;
		demandsum_crop=0.0;
		supplysum_crop=0.0;

		growingseason=false;	//Initialized to true for normal grass growth (CC3G & CC4G) in establishment
		growingseason_ystd=false;
		senescence=false;
		senescence_ystd=false;
		intercropseason=false;
		bicdate=-1;
		eicdate=-1;
		growingdays=0;
		growingdays_y=0;
		lgp=0;

		for(int j=0;j<2;j++) {
			sdate_harvest[j]=-1;
			hdate_harvest[j]=-1;
			sdate_thisyear[j]=-1;
		}

		vdsum_alloc=0.0;
		vd = 0.0;
		f_alloc_root=0.0;
		f_alloc_leaf=0.0;
		f_alloc_horg=0.0;
		f_alloc_stem=0.0;
		dev_stage = 0.0;

		fertilised[0] = false;
		fertilised[1] = false;
		fertilised[2] = false;
	}

	void serialize(ArchiveStream& arch);
};


/// State variables common to all individuals of a particular PFT in a particular patch
/** Used in individual and cohort modes only. */
class Patchpft : public Serializable {

	// MEMBER VARIABLES:

public:

	/// id code (equal to value of member variable id in corresponding Pft object)
	int id;
	/// reference to corresponding Pft object in PFT list
	Pft& pft;
	/// potential annual net assimilation (leaf-level net photosynthesis) at forest floor (kgC/m2/year)
	double anetps_ff;
	/// water stress parameter (0-1 range; 1=minimum stress)
	double wscal;
	/// running sum (converted to annual mean) for wscal
	double wscal_mean;
	/// potential annual net assimilation at forest floor averaged over establishment interval (kgC/m2/year)
	double anetps_ff_est;
	/// first-year value of anetps_ff_est
	double anetps_ff_est_initial;
	/// annual mean wscal averaged over establishment interval
	double wscal_mean_est;
	/// vegetation phenological state (fraction of potential leaf cover), updated daily
	double phen;
	/// monthly sum of daily fractional leaf cover
	double mphen[12];
	/// annual sum of daily fractional leaf cover
	/** equivalent number of days with full leaf cover
	 *  (reset on expected coldest day of year)
	 */
	double aphen;
	/// Month of lowest phen for RAINGREEN
	int driest_mth;
	/// whether PFT can establish in this patch under current conditions
	bool establish;
	/// running total for number of saplings of this PFT to establish (cohort mode)
	double nsapling;
	/// leaf-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_leaf;
	/// fine root-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_root;
	/// remaining sapwood-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_sap;
	/// remaining heartwood-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_heart;
	/// litter derived from allocation to reproduction for PFT on modelled area basis (kgC/m2)
	double cmass_litter_repr;

	/// leaf-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_leaf;
	/// root-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_root;
	/// remaining sapwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_sap;
	/// remaining heartwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_heart;

	/// leaf-derived phosphorus litter for PFT on modelled area basis (kgP/m2)
	double pmass_litter_leaf;
	/// root-derived phosphorus litter for PFT on modelled area basis (kgP/m2)
	double pmass_litter_root;
	/// remaining sapwood-derived phosphorus litter for PFT on modelled area basis (kgP/m2)
	double pmass_litter_sap;
	/// remaining heartwood-derived phosphorus litter for PFT on modelled area basis (kgP/m2)
	double pmass_litter_heart;

	/// non-FPC-weighted canopy conductance value for PFT under water-stress conditions (mm/s)
	double gcbase;
	/// daily value of the above variable (mm/s)
	double gcbase_day;

	/// evapotranspirational "supply" function for this PFT today (mm/day)
	double wsupply;
	double wsupply_leafon;
	/// fractional uptake of water from each soil layer today
	double fwuptake[NSOILLAYER];

	/// whether water-stress conditions for this PFT
	bool wstress;
	/// daily version of the above variable
	bool wstress_day;

	/// carbon depository for long-lived products like wood
	double cmass_harvested_products_slow;
	/// nitrogen depository for long-lived products like wood
	double nmass_harvested_products_slow;
	/// phosphorus depository for long-lived products like wood
	double pmass_harvested_products_slow;

	/// first and last day of crop sowing window, calculated in crop_sowing_patch() or Crop_sowing_date_new()
	int swindow[2];
	/// daily value of water deficit, calculated in irrigated_water_uptake()
	double water_deficit_d;
	/// yearly sum of water deficit
	double water_deficit_y;

	/// Struct for crop-specific variables
	cropphen_struct *cropphen;

	// INUNDATION STRESS TERMS:
    
	/// Number of days a month that the water table is above this PFT's wtp_max (updated on the first day of the month)
	int inund_count;
	/// [0,1] - a measure of the inundation stress. Daily photosynthesis is reduced by this factor.
	double inund_stress;


	// MEMBER FUNCTIONS:

	/// Constructor: initialises id, pft and data members
	Patchpft(int i,Pft& p):id(i),pft(p) {

		cmass_litter_leaf = 0.0;
		cmass_litter_root = 0.0;
		cmass_litter_sap   = 0.0;
		cmass_litter_heart = 0.0;
		cmass_litter_repr = 0.0;

		nmass_litter_leaf  = 0.0;
		nmass_litter_root  = 0.0;
		nmass_litter_sap   = 0.0;
		nmass_litter_heart = 0.0;

		pmass_litter_leaf = 0.0;
		pmass_litter_root = 0.0;
		pmass_litter_sap = 0.0;
		pmass_litter_heart = 0.0;

		wscal = 1.0;
		wscal_mean = 1.0;
		anetps_ff = 0.0;
		aphen = 0.0;
		phen = 0.0;
		driest_mth = 0;
		wsupply = 0.0;
		wsupply_leafon = 0.0;
		anetps_ff_est = 0.0;
		anetps_ff_est_initial = 0.0;
		wscal_mean_est = 0.0;
		nsapling = 0;

		for (int mth = 0; mth < 12; mth++)
			mphen[mth] = 0.0;

		for (int i = 0; i < NSOILLAYER; i++)
			fwuptake[i] = 0.0;

		cropphen = NULL;
		cmass_harvested_products_slow = 0.0;
		nmass_harvested_products_slow = 0.0;
		pmass_harvested_products_slow = 0.0;

		swindow[0]=-1;
		swindow[1]=-1;

		if(pft.landcover==CROPLAND)
		{
			cropphen=new cropphen_struct;
		}

		inund_count=0;
		inund_stress=1.0; // No stress by default
	}

	~Patchpft() {
		if(cropphen) {
			delete cropphen;
		}
	}

	/// safe method to obtain cropphen_struct pointer
	cropphen_struct* get_cropphen();
	/// safe method to set cropphen_struct variables
	cropphen_struct* set_cropphen();

	bool growingseason() const;

	void serialize(ArchiveStream& arch);
};


/// Stores data for a patch.
/** In cohort and individual modes, replicate patches are
 *  required in each stand to accomodate stochastic variation; in population mode there
 *  should be just one Patch object, representing average conditions for the entire
 *  stand. A reference to the parent Stand object (defined below) is included as a
 *  member variable.
 */
class Patch : public Serializable {

public:

	// MEMBER VARIABLES

	/// id code in range 0-npatch for patch
	int id;
	/// reference to parent Stand object
	Stand& stand;
	/// list array [0...npft-1] of Patchpft objects (initialised in constructor)
	ListArray_idin1<Patchpft,Pft> pft;
	/// vegetation for this patch
	Vegetation vegetation;
	/// soil for this patch
	Soil soil;
	/// fluxes for this patch
	Fluxes fluxes;
	/// FPAR at top of grass canopy today
	double fpar_grass;
	/// FPAR at soil surface today
	double fpar_ff;
	/// mean growing season PAR at top of grass canopy (J/m2/day)
	double par_grass_mean;
	/// number of days in growing season, estimated from mean vegetation leaf-on fraction
	/** \see function fpar in canopy exchange module */
	int nday_growingseason;
	/// total patch FPC
	double fpc_total;
	/// whether patch was disturbed last year
	bool disturbed;
	/// patch age (years since last disturbance)
	int age;
	/// probability of fire this year (GlobFIRM)
	double fireprob;

	/// BLAZE if patch has burnt this year
	bool burned;
	/// BLAZE Fire line intensity (kW/m)
	double fire_line_intensity;

	// BLAZE fire related carbon fluxes
	/// BLAZE-fire carbon flux: live wood to atmosphere (kgC/m2)
	double wood_to_atm;
	/// BLAZE-fire carbon flux: leaves to atmosphere (kgC/m2)
	double leaf_to_atm;
	/// BLAZE-fire carbon flux: leaves to litter (kgC/m2)
	double leaf_to_lit;
	/// BLAZE-fire carbon flux: live wood to structural litter (kgC/m2)
	double wood_to_str;
	/// BLAZE-fire carbon flux: live wood to fine woody debris (kgC/m2)
	double wood_to_fwd;
	/// BLAZE-fire carbon flux: live wood to coarse woody debris (kgC/m2)
	double wood_to_cwd;
	/// BLAZE-fire carbon flux: fine litter to atmosphere (kgC/m2)
	double litf_to_atm;
	/// BLAZE-fire carbon flux: fine woody debris to atmosphere (kgC/m2)
	double lfwd_to_atm;
	/// BLAZE-fire carbon flux: coarse woody debris to atmosphere (kgC/m2)
	double lcwd_to_atm;

	// Storage for averaging of different Fapars for biome mapping in SIMFIRE
	/// SIMFIRE fapar: Grasses
	double fapar_grass_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Needle-leaf tree
	double fapar_ndlt_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Broad-leaf tree
	double fapar_brlt_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: TrBR (Savanna)
	double fapar_trbr_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Shrubs
	double fapar_shrub_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Total fapar
	double fapar_total_avg[N_YEAR_BIOMEAVG];

	/// whether management has started on this patch
	bool managed;
	/// cutting intensity (initial percent of trees cut, further selection at individual level has to be done in a separate function)
	double man_strength;

	bool managed_this_year;
	bool plant_this_year;

	/// DLE - the number of days over which wcont is averaged for this patch
	/** i.e. those days for which daily temp > 5.0 degC */
	int growingseasondays;


	// Variables used by new hydrology (Dieter Gerten 2002-07)

	/// interception by vegetation today on patch basis (mm)
	double intercep;
	/// annual sum of AET (mm/year)
	double aaet;
	/// annual sum of AET (mm/year) for each of the last five simulation years
	Historic<double, NYEARAAET> aaet_5;
	/// annual sum of soil evaporation (mm/year)
	double aevap;
	/// annual sum of interception (mm/year)
	double aintercep;
	/// annual sum of runoff (mm/year)
	double asurfrunoff;
	/// annual sum of runoff (mm/year)
	double adrainrunoff;
	/// annual sum of runoff (mm/year)
	double abaserunoff;
	/// annual sum of runoff (mm/year)
	double arunoff;
	/// water added to wetlands today (mm)
	double wetland_water_added_today;
	/// annual sum of water added to wetlands (mm/year)
	double awetland_water_added;
	/// annual sum of potential evapotranspiration (mm/year)
	double apet;

	/// equilibrium evapotranspiration today, deducting interception (mm)
	double eet_net_veg;

	/// transpirative demand for patch, patch vegetative area basis (mm/day)
	double wdemand;
	/// daily average of the above variable (mm/day)
	double wdemand_day;
	/// transpirative demand for patch assuming full leaf cover today
	/** mm/day, patch vegetative area basis	*/
	double wdemand_leafon;
	/// rescaling factor to account for spatial overlap between individuals/cohorts populations
	double fpc_rescale;

	/// monthly AET (mm/month)
	double maet[12];
	/// monthly soil evaporation (mm/month)
	double mevap[12];
	/// monthly interception (mm/month)
	double mintercep[12];
	/// monthly runoff (mm/month)
	double mrunoff[12];
	/// monthly PET (mm/month)
	double mpet[12];

	/// daily nitrogen demand
	double ndemand;

	/// daily phosphorus demand
	double pdemand;

	/// annual nitrogen fertilization (kgN/m2/year)
	double anfert;
	/// daily nitrogen fertilization (kgN/m2/day)
	double dnfert;

	/// annual phosphorus fertilization (kgP/m2/year)
	double apfert;
	/// daily phosphorus fertilization (kgP/m2/day)
	double dpfert;

	/// daily value of irrigation water (mm), set in irrigation(), derived from water_deficit_d
	double irrigation_d;
	/// yearly sum of irrigation water (mm)
	double irrigation_y;

	/// whether litter is to be sent to the soil today
	bool is_litter_day;
	/// number of harvests and/or cover-crop killing or turnover events
	int nharv;
	/// whether today is a harvest day and/or cover-crop killing or turnover day
	bool isharvestday;

	// MEMBER FUNCTIONS

	/// Constructor: initialises various members and builds list array of Patchpft objects.
	Patch(int i,Stand& s,Soiltype& st);

	void serialize(ArchiveStream& arch);

	/// Returns the Climate for this Patch
	/** This function returns a const reference to prevent code which operates
	 *  on a patch basis to modify the climate and thereby affect other
	 *  patches/stands.
	 */
	const Climate& get_climate() const;

	/// Returns whether we should model fire in this patch
	bool has_fires() const;

	/// Returns whether we should model disturbances in this patch
	bool has_disturbances() const;
	/// Total patch carbon biomass and litter
	double ccont(double scale_indiv = 1.0, bool luc = false);
	/// Total patch nitrogen biomass and litter
	double ncont(double scale_indiv = 1.0, bool luc = false);
	/// Total patch phosphorus biomass and litter
	double pcont(double scale_indiv = 1.0, bool luc = false);
	/// Total patch carbon fluxes so far this year
	double cflux();
	/// Total patch nitrogen fluxes so far this year
	double nflux();
	/// Total patch phosphorus fluxes so far this year
	double pflux();
	
	/// Get 5-year mean of wood C mass increase (periodic annual increment)
	double get_cmass_wood_inc_5() {
		double cmass_wood_inc_5_mean = 0.0;
		for (unsigned int i=0; i<vegetation.nobj; i++) {

			// Disregard shrubs (crownarea_max = 10)
			Individual& indiv = vegetation[i];
			if(indiv.pft.lifeform == TREE && indiv.pft.crownarea_max > 10) {
				if(indiv.cmass_wood_inc_5.size())
					cmass_wood_inc_5_mean += indiv.cmass_wood_inc_5.mean();
			}
		}
		return cmass_wood_inc_5_mean;
	}
	
	/// Get cmass_wood of all individuals in patch
	double cmass_wood() {
		double cmass_wood = 0.0;
		for (unsigned int i=0; i<vegetation.nobj; i++) {

			Individual& indiv = vegetation[i];
			cmass_wood += indiv.cmass_wood();
		}
		return cmass_wood;
	}
};

/// Container for variables common to individuals of a particular PFT in a stand.
/** Used in individual and cohort modes only
 */
class Standpft : public Serializable {

public:

	// MEMBER VARIABLES

	int id;
	Pft& pft;
	/// net C allocated to reproduction for this PFT in all patches of this stand this year (kgC/m2)
	double cmass_repr;
	/// maximum value of Patchpft::anetps_ff for this PFT in this stand so far in the simulation (kgC/m2/year)
	double anetps_ff_max;

	/// FPC sum for this PFT as average for stand
	double fpc_total;

	/// Photosynthesis values for this PFT under non-water-stress conditions
	PhotosynthesisResult photosynthesis;

	/// Whether this PFT is allowed to grow in this stand
	bool active;
	/// Whether this PFT is planted in this stand
	bool plant;
	/// Whether this PFT is allowed to establish (after planting) in this stand
	bool reestab;

	/// Whether this PFT is irrigated in this stand
	bool irrigated;
	/// sowing date specified in stand type or read from input file
	int sdate_force;
	/// harvest date specified in stand type or read from input file
	int hdate_force;

	// MEMBER FUNCTIONS

	/// Constructor: initialises various data members
	Standpft(int i,Pft& p):id(i),pft(p) {

		anetps_ff_max = 0.0;
		active = !run_landcover;
		plant = false;
		reestab = false;
		irrigated = false;
		sdate_force = -1;
		hdate_force = -1;
	}

	void serialize(ArchiveStream& arch);
};


/// The stand class corresponds to a modelled area of a specific landcover type in a grid cell.
/** There may be several stands of the same landcover type (but with different settings).
 */
class Stand : public ListArray_idin2<Patch,Stand,Soiltype>, public Serializable {

public:

	// MEMBER VARIABLES

	/// list array [0...npft-1] of Standpft (initialised in constructor)
	ListArray_idin1<Standpft,Pft> pft;

	/// A number identifying this Stand within the grid cell
	int id;

	/// stand type id
	int stid;

	/// pft id of main crop, updated during rotation
	int pftid;

	/// current crop rotation item
	int current_rot;
	/// number of days passed in current rotation item
	int ndays_inrotation;
	/// Returns true if stand is in fallow (with cover crop grass)
	bool infallow;
	/// Returns true if crop rotation item is to be updated today
	bool isrotationday;
	/// Returns true if current crop management hydrology == irrigated, updated during rotation
	bool isirrigated;
	/// Returns true if the stand's main crop pft intercrop==naturalgrass and a pft with isintercrop==true is in the pftlist.
	bool hasgrassintercrop;
	/// gdd5-value at first intercrop grass growth
	double gdd5_intercrop;

	/// old fraction of this stand relative to the gridcell before update
	double frac_old;

	/// used during land cover change involving several calls to reveiving_stand_change()
	/** Set to frac_old in reduce_stands(), then modified in donor_stand_change() and receiving_stand_change().
	 */
	double frac_temp;
	/// fraction unavailable for transfer to other stand types
	double protected_frac;
	/// net stand fraction change
	double frac_change;
	/// gross fraction increase
	double gross_frac_increase;
	/// gross fraction decrease
	double gross_frac_decrease;
	/// fraction that has been cloned from another stand
	double cloned_fraction;
	/// Returns true if this stand is cloned from another stand
	bool cloned;
	/// pointer to array of fractions transferred from this stand to other stand types
	double *transfer_area_st;
	/// land cover origin of this stand
	landcovertype origin;
	/// used for output from separate stands
	double anpp;
	/// used for output from separate stands
	double cmass;

	/// Seed for generating random numbers within this Stand
	/** The reason why Stand has its own seed, rather than using for instance
	 *  a single global seed is to make it easier to compare results when using
	 *  different land cover types.
	 *
	 *  Randomness not associated with a specific stand, but rather a whole
	 *  grid cell should instead use the seed in the Gridcell class.
	 *
	 *  \see randfrac()
	 */
	long seed;

	/// type of landcover
	/** \see landcovertype
	 *  initialised in constructor
	 */
	landcovertype landcover;

	/// The year when this stand was created.
	/** Will typically be year zero unless running with dynamic
	 *  land cover.
	 *
	 *  Needed to set patchpft.anetps_ff_est_initial
	 */
	int first_year;
	// The year this stand was cloned from another stand
	int clone_year;
	/// scaling factor for stands that have grown in area this year (old fraction/new fraction)
	double scale_LC_change;

	// MEMBER FUNCTIONS

	/// Constructs a Stand
	/** \param i         The id for the stand within the grid cell
	 *  \param gc        The parent grid cell
	 *  \param st        The soil type to be used within this Stand
	 *  \param landcover The type of landcover to use for this stand
	 */
	Stand(int i, Gridcell* gc, Soiltype& st, landcovertype landcover, int no_patch = 0);

	~Stand();

	/// Gives the fraction of this Stand relative to the whole grid cell
	double get_gridcell_fraction() const;

	/// Gives the fraction of this Stand relative to its land cover type; NB: unsafe to use within landcover_dynamics() !
	double get_landcover_fraction() const;

	/// Set the fraction of this Stand relative to the gridcell
	void set_gridcell_fraction(double fraction);

	/// Returns the number of patches in this Stand
	unsigned int npatch() const { return nobj; }

	/// Returns true if stand is true crop stand, as opposed to pasture grass grown on cropland or other land cover
	inline bool is_true_crop_stand() {
		return landcover==CROPLAND && pft[pftid].pft.phenology==CROPGREEN;	// OK also for fallow (pftid always cropgreen)
	}
	/// Moves crop rotation forward
	void rotate();
	/// Returns area transferred to other land cover during land cover change
	double transfer_area_lc(landcovertype to);
	/// Initiates new stand land cover settings
	void init_stand_lu(StandType& st, double fraction);
	/// Total stand carbon biomass and litter
	double ccont(double scale_indiv = 1.0);
	/// Total stand nitrogen biomass and litter
	double ncont(double scale_indiv = 1.0);
	/// Total stand phosphorus biomass and litter
	double pcont(double scale_indiv = 1.0);
	/// Total stand carbon fluxes so far this year
	double cflux();
	/// Total stand nitrogen fluxes so far this year
	double nflux();
	/// Total stand phosphorus fluxes so far this year
	double pflux();
	/// Returns true if stand is true high-latitude peatland stand, as opposed to a wetland < PEATLAND_WETLAND_LATITUDE_LIMIT N
	bool is_highlatitude_peatland_stand() const;
	/// Returns true if stand is wetland stand, as opposed to a peatland >= PEATLAND_WETLAND_LATITUDE_LIMIT N
	bool is_true_wetland_stand() const;
	
	/// Creates a duplicate stand with a new landcovertype
	/** The new stand is added to this stand's gridcell.
	*
	*  \returns reference to the new stand
	*/
	Stand& clone(StandType& st, double fraction);

	void serialize(ArchiveStream& arch);

	/// Returns the Climate for this Stand
	/** This function returns a const reference to prevent code which operates
	 *  on a stand basis to modify the climate and thereby affect other
	 *  stands.
	 */
	const Climate& get_climate() const;

	/// Returns the Gridcell containing this Stand
	Gridcell& get_gridcell() const;

private:

	/// Pointer to parent object, could be a null pointer
	/** Prefer to access the gridcell through get_gridcell(), even internally
	 *  within the Stand class.
	 */
	Gridcell* gridcell;

	/// Soil type to be used in this Stand
	Soiltype& soiltype;

	/// Fraction of this stand relative to the gridcell
	/** used by crop stands; initialized in constructor to 1,
	 *  set in landcover_init()
	 */
	double frac;
};



/// State variables common to all individuals of a particular PFT in a GRIDCELL.
class Gridcellpft : public Serializable {

public:

	// MEMBER VARIABLES

	/// A number identifying this object within its list array
	int id;

	/// A reference to the Pft object for this Gridcellpft
	Pft& pft;

	/// annual degree day sum above threshold damaging temperature
	/** used in calculation of heat stess mortality; Sitch et al 2000, Eqn 55
	 */
	double addtw;

	/// Michaelis-Menten kinetic parameters
	/** Half saturation concentration for N uptake (Rothstein 2000, Macduff 2002)
	 */
	double Km;
	/** Half saturation concentration for P uptake (Silveira & Cardoso 2004)
	*/
	double Kmp;

	///Crop-specific variables:
	/// whether the daily temperature has fallen below the autumn temperature limit (tempautumn) this year
	bool autumnoccurred;
	/// whether the daily temperature has risen above the spring temperature limit (tempspring) this year
	bool springoccurred;
	/// whether the daily temperature has fallen below the vernalization limit (trg) this year
	bool vernstartoccurred;
	/// whether the daily temperature rises over the vernalization limit (trg) this year
	bool vernendoccurred;
	/// first day when temperature fell below the autumn temperature limit (tempautumn) this year
	int first_autumndate;
	/// 20-year mean
	int first_autumndate20;
	/// memory of the last 20 years' values
	int first_autumndate_20[20];
	/// last day when temperature rose above the spring temperature limit (tempspring) this year
	int last_springdate;
	/// 20-year mean
	int last_springdate20;
	/// memory of the last 20 years' values
	int last_springdate_20[20];
	/// last day when temperature has fallen below the vernilisation temperature limit (trg) this year (if vernstartoccurred==true)
	int last_verndate;
	/// 20-year mean
	int last_verndate20;
	/// memory of the last 20 years' values
	int last_verndate_20[20];
	/// default sowing date (pft.sdatenh/sdatesh)
	int sdate_default;
	/// calculated sowing date from temperature limits
	int sdatecalc_temp;
	/// calculated sowing date from precipitation limits
	int sdatecalc_prec;
	/// sowing date from input file
	int sdate_force;
	/// harvest date from input file
	int hdate_force;
	/// N fertilization from input file
	double Nfert_read;
	/// Manure N fertilization from input file
	double Nfert_man_read;
	/// P fertilization from input file
	double Pfert_read;
	/// Manure P fertilization from input file
	double Pfert_man_read;
	/// default harvest date (pft.hlimitdatenh/hlimitdatesh)
	int hlimitdate_default;
	/// whether autumn sowing is either calculated or prescribed
	bool wintertype;
	/// first and last day of crop sowing window, calculated in calc_sowing_windows()
	int swindow[2];
	/// first and last day of crop sowing window for irrigated crops, calculated in calc_sowing_windows()
	int swindow_irr[2];
	/// temperature limits precludes crop sowing
	bool sowing_restriction;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcellpft object
	/** \param i   The id for this object
	 *  \param p   A reference to the Pft for this Gridcellpft
	 */
	Gridcellpft(int i,Pft& p):id(i),pft(p) {
		addtw = 0.0;
		Km = 0.0;

		autumnoccurred=false;
		springoccurred=false;
		vernstartoccurred=false;
		vernendoccurred=false;
		first_autumndate=-1;
		first_autumndate20=-1;
		last_springdate=-1;
		last_springdate20=-1;
		last_verndate=-1;
		last_verndate20=-1;
		for (int year=0;year<20;year++) {
			first_autumndate_20[year]=-1;
			last_springdate_20[year]=-1;
			last_verndate_20[year]=-1;
		}
		sdate_default=-1;
		sdate_force=-1;
		hdate_force=-1;
		Nfert_read=-1;
		Nfert_man_read=-1;
		Pfert_read = -1;
		Pfert_man_read = -1;
		sdatecalc_temp=-1;
		sdatecalc_prec=-1;
		hlimitdate_default=-1;
		wintertype=false;
		swindow[0]=-1;
		swindow[1]=-1;
		sowing_restriction = false;
	}

	void serialize(ArchiveStream& arch);
};

/// State variables common to all individuals of a particular STANDTYPE in a GRIDCELL.
class Gridcellst : public Serializable {

public:

	// MEMBER VARIABLES

	/// A number identifying this object within its list array
	int id;

	/// A reference to the StandType object for this Gridcellst
	StandType& st;

	/// fraction of this stand type relative to the gridcell
	double frac;
	/// old fraction of this stand type relative to the gridcell before update
	double frac_old;
	/// original input value of old fraction of this stand type before rescaling
	double frac_old_orig;
	/// fraction unavailable for transfer to other stand types
	double protected_frac;

	/// net fraction change
	double frac_change;
	/// gross fraction increase
	double gross_frac_increase;
	/// gross fraction decrease
	double gross_frac_decrease;

	// current number of stands of this stand type
	int nstands;

	double nfert;

	double pfert;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcellst object
	/** \param i   The id for this object
	 *  \param s   A reference to the StandType for this Gridcellst
	 */
	Gridcellst(int i,StandType& s):id(i),st(s) {
		frac = 1.0;
		frac_old = 0.0;
		frac_old_orig = 0.0;
		protected_frac = 0.0;
		frac_change = 0.0;
		gross_frac_increase = 0.0;
		gross_frac_decrease = 0.0;
		nstands = 0;
		nfert = -1.0;
		pfert = -1.0;
	}

	void serialize(ArchiveStream& arch);
};

/// Storage of land cover fraction data and some land cover change-related pools and fluxes
struct Landcover : public Serializable {

	Landcover();

	/// The fractions of the different land cover types.
	/** landcoverfrac is read in from land cover input file or from
	 *  instruction file in getlandcover().
	 */
	double frac[NLANDCOVERTYPES];

	/// The land cover fractions from the previous year
	/** Used to keep track of the changes when running with dynamic
	 *  land cover.
	 */
	double frac_old[NLANDCOVERTYPES];

	double frac_change[NLANDCOVERTYPES];

	/// Transfer matrices
	double frac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES];
	double primary_frac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES];

	/// Whether the land cover fractions changed for this grid cell this year
	/** \see landcover_dynamics
	 */
	bool updated;

	/// Gridcell-level C flux from slow harvested products
	double acflux_harvest_slow;

	/// Gridcell-level C flux from harvest associated with landcover change
	double acflux_landuse_change;

	/// Gridcell-level N flux from slow harvested products
	double anflux_harvest_slow;

	/// Gridcell-level N flux from harvest associated with landcover change
	double anflux_landuse_change;

	/// Gridcell-level P flux from slow harvested products
	double apflux_harvest_slow;

	/// Gridcell-level P flux from harvest associated with landcover change
	double apflux_landuse_change;

	/// Landcover-level C flux from slow harvested products (donating landcover)
	double acflux_harvest_slow_lc[NLANDCOVERTYPES];

	/// Landcover-level C flux from harvest associated with landcover change (donating landcover)
	double acflux_landuse_change_lc[NLANDCOVERTYPES];

	/// Landcover-level N flux from slow harvested products (donating landcover)
	double anflux_harvest_slow_lc[NLANDCOVERTYPES];

	/// Landcover-level N flux from harvest associated with landcover change (donating landcover)
	double anflux_landuse_change_lc[NLANDCOVERTYPES];

	/// Landcover-level P flux from slow harvested products (donating landcover)
	double apflux_harvest_slow_lc[NLANDCOVERTYPES];

	/// Landcover-level P flux from harvest associated with landcover change (donating landcover)
	double apflux_landuse_change_lc[NLANDCOVERTYPES];

	/// Which landcover types create new stands when area increases.
	bool expand_to_new_stand[NLANDCOVERTYPES];

	/// Whether to pool all transferred land from a donor landcover (overrides different landcover targets of different stand types and stands in a landcover)
	bool pool_to_all_landcovers[NLANDCOVERTYPES];

	/// Whether to pool transferred land to a receptor landcover (crop and pasture stands to new natural stand: pool!)
	bool pool_from_all_landcovers[NLANDCOVERTYPES];

	void serialize(ArchiveStream& arch);
};

/// The Gridcell class corresponds to a modelled locality or grid cell.
/** Member variables include an object of type Climate (holding climate, insolation and
 *  CO2 data), a object of type Soiltype (holding soil static parameters) and a list
 *  array of Stand objects. Soil objects (holding soil state variables) are associated
 *  with patches, not gridcells. A separate Gridcell object must be declared for each modelled
 *  locality or grid cell.
 */
class Gridcell : public GuessContainer<Stand>, public Serializable {

public:

	// MEMBER VARIABLES

	/// climate, insolation and CO2 for this grid cell
	Climate climate;

	/// soil static parameters for this grid cell
	Soiltype soiltype;

	/// landcover fractions and landcover-specific variables
	Landcover landcover;

	/// list array [0...npft-1] of Gridcellpft (initialised in constructor)
	ListArray_idin1<Gridcellpft,Pft> pft;

	/// list array [0...nst-1] of Gridcellst (initialised in constructor)
	ListArray_idin1<Gridcellst,StandType> st;

	/// object for keeping track of carbon and nitrogen balance
	MassBalance balance;

	// SIMFIRE
	/// the region index to chosose from set of optimisations
	int simfire_region;
	/// timeseries of population density from the Hyde 3.1 dataset (inhabitants/ha)
	double hyde31_pop_density[57];
	/// current year's population density (inhabitants/ha)
	double pop_density;
	/// tuning factor for available litter
	double k_tun_litter;
	/// maximum annual Nesterov Index
	double nesterov_max;
	/// current Nexterov index
	double nesterov_cur;
	/// Monthly max Nexterov index (to keep track of running year)
	double nesterov_monthly_max[12];
	/// biome classification used in SIMFIRE
	int simfire_biome;
	/// Average maximum annual fAPAR (over avg_interv_fpar years)
	double fapar_annual_max;
	/// Average maximum annual fAPAR of recent years
	double fapar_recent_max[AVG_INTERVAL_FAPAR];
	/// maximum fapar of running year so far
	double fapar_cur_max;
	/// monthly fire risk (factor describing local monthly fire climatology)
	double monthly_fire_risk[12];
	/// current burned area from SIMFIRE (fract.)
	double burned_area;
	/// effectively burned area (after randFrac in BLAZE)
	double effective_burned_area;
	/// accumulated burned area from SIMFIRE for tstep < 1a (fract.)
	double simfire_annual_burned_area;
	/// Simple tracker to check whether at least one patch has enough fuel to burn
	int can_burn;
	/// annual burned area
	double annual_burned_area;
	/// monthly burned area
	double monthly_burned_area[12];

	// Nitrogen deposition
	/// annual NH4 deposition (kgN/m2/year)
	double aNH4dep;
	/// annual NO3 deposition (kgN/m2/year)
	double aNO3dep;
	/// annual P deposition (kgP/m2/year)
	double apdep;
	/// daily NH4 deposition (kgN/m2)
	double dNH4dep;
	/// daily NO3 deposition (kgN/m2)
	double dNO3dep;
	/// daily P deposition (kgP/m2)
	double dpdep;

	/// Seed for generating random numbers within this Gridcell
	/** The reason why Gridcell has its own seed, rather than using for instance
	 *  a single global seed is to make it easier to compare results when for
	 *  instance changing the order in which the simulation proceeds. It also
	 *  gets serialized together with the rest of the Gridcell state to make it
	 *  possible to get exactly identical results after a restart.
	 *
	 *  \see randfrac()
	 */
	long seed;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcell object
	Gridcell();

	/// Longitude for this grid cell
	double get_lon() const;

	/// Latitude for this grid cell
	double get_lat() const;

	/// Set longitude and latitude for this grid cell
	void set_coordinates(double longitude, double latitude);

	void serialize(ArchiveStream& arch);

	/// Creates a new Stand in this grid cell
	Stand& create_stand(landcovertype lc, int no_patch = 0);

	/// Creates new stand and initiates land cover settings when run_landcover==true
	Stand& create_stand_lu(StandType& st, double fraction, int no_patch = 0);

	/// Total gridcell carbon biomass and litter
	double ccont();
	/// Total gridcell nitrogen biomass and litter
	double ncont();
	/// Total gridcell phosphorus biomass and litter
	double pcont();
	/// Total gridcell carbon fluxes so far this year
	double cflux();
	/// Total gridcell nitrogen fluxes so far this year
	double nflux();
	/// Total gridcell phosphorus fluxes so far this year
	double pflux();

	/// Deletes the stand which the iterator is pointing at
	/** Returns an iterator pointing to the object following the erased object.
	 */
	iterator delete_stand(iterator itr);

	/// Returns number of stands
	unsigned int nbr_stands() const;

private:

	/// Longitude for this grid cell
	double lon;

	/// Latitude for this grid cell
	double lat;

};


#endif // LPJ_GUESS_GUESS_H

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Delmas, R., Lacaux, J.P., Menaut, J.C., Abbadie, L., Le Roux, X., Helaa, G., Lobert, J., 1995.
//   Nitrogen compound emission from biomass burning in tropical African Savanna FOS/DECAFE 1991
//   experiment. Journal of Atmospheric Chemistry 22, 175-193.
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical
//   exploration of the relationships of soil moisture characteristic to the
//   physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Franzlubbers, AJ & Stuedemann, JA 2009 Soil-profile organic carbon and total
//   nitrogen during 12 years of pasture management in the Southern Piedmont USA.
//   Agriculture Ecosystems & Environment, 129, 28-36.
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A
//   process-based, terrestrial biosphere model of ecosystem dynamics
//   (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Fulton, MR 1991 Adult recruitment rate as a function of juvenile growth in size-
//   structured plant populations. Oikos 61: 102-105.
// Gerten, D., Schaphoff, S., Haberlandt, W., Lucht, W. & Sitch, S. 2004. 
//   Terrestrial vegetation and water balance—hydrological evaluation of a dynamic 
//   global vegetation model. Journal of Hydrology 286: 249-270.
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Jackson, R.B., Canadell, J., Ehleringer, J.R., Mooney, H.A., Sala O.E. & Schulze, E.D. 1996
//   A global analysis of root distributions for terrestrial biomes. 
//   Oecologia, Volume 108: 389–411
// Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
//   Functional Ecology 8: 315-323
// Macduff, JH, Humphreys, MO & Thomas, H 2002. Effects of a stay-green mutation on
//   plant nitrogen relations in Lolium perenne during N starvation and after
//   defoliation. Annals of Botany, 89, 11-21.
// Monsi M & Saeki T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften und
//   seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany 14: 22-52
// Olin S., G. Schurgers, M. Lindeskog, D. Wårlind, B. Smith, P. Bodin, J.
//    Holmér, and A. Arneth. 2015 Biogeosciences Discuss., 12, 1047-1111. The
//    impact of atmospheric CO2 and N management on yields and tissue C:N in
//    the main wheat regions of Western Europe
// Parton, W. J., Hanson, P. J., Swanston, C., Torn, M., Trumbore, S. E., Riley, W.
//   & Kelly, R. 2010. ForCent model development and testing using the Enriched
//   Background Isotope Study experiment. Journal of Geophysical
//   Research-Biogeosciences, 115.
// Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65: 51-70.
// Reich, PB, Walters MB & Ellsworth DS 1992 Leaf Life-Span in Relation to Leaf,
//   Plant, and Stand Characteristics among Diverse Ecosystems.
//   Ecological Monographs 62: 365-392.
// Silveira, A. P. D., Cardoso, E. J. B. N. (2004) Arbuscular mycorrhiza and kinetic 
//	 parameters of phosphorus absorption by bean plants.Sci.Agric. 
//	 (Piracicaba, Braz.), v.61, n.2, p.203 - 209
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Sykes, MT, Prentice IC & Cramer W 1996 A bioclimatic model for the potential
//   distributions of north European tree species under present and future climates.
//   Journal of Biogeography 23: 209-233.
// Wang, E, Engel, T, 1998 Simulation of phenological development of wheat crops
//   Agricultural Systems 58:1-24
// White, M A, Thornton, P E, Running, S. & Nemani, R 2000 Parameterization and
//   Sensitivity Analysis of the BIOME-BGC Terrestrial Ecosystem Model: Net Primary
//   Production Controls. Earth Interactions, 4, 1-55.
