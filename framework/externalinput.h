////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file externalinput.h
/// \brief Input code for land cover, management and other data, currently from text files.
/// \author Mats Lindeskog
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_EXTERNALINPUT_H
#define LPJ_GUESS_EXTERNALINPUT_H

/// Landcover area fraction input resolution used in the code to reject changes caused by rounding errors.
extern double INPUT_RESOLUTION;

#include "indata.h"

using namespace TextInput;

/// Reads gridlist in lon-lat-description format from text input file
void read_gridlist(ListArray_id<Coord>& gridlist, const char* file_gridlist);

/// Help function for get_lc_transfer() to adjust inconsistencies between net land cover inout and gross land cover transitions.
void adjust_gross_transfers(Gridcell& gridcell, double landcoverfrac_change[], double lc_frac_transfer[][NLANDCOVERTYPES],
							forest_lc_frac_transfer& forest_lc_subset_transfer, double& tot_frac_change);

/// Class that deals with miscellaneous data input.
/** Input data that are not from the input module or from the land-cover or management input modules.
 */
class MiscInput {

public:

	/// Constructor
	MiscInput() {;}

	/// Opens land cover input files
	void init();

	/// Loads disturbance from input file
	bool loaddisturbance(double lon, double lat);

	/// Loads elevation input file
	bool loadelevation(double lon, double lat);

	/// Gets all static input data from this class
	void getmiscinput_static(Gridcell& gridcell);

	/// Gets all yearly input data from this class
	void getmiscinput_yearly(Gridcell& gridcell);

private:

	// Objects handling additional environmental data input
	TextInput::TimeDataD disturbance;
	TextInput::TimeDataD disturbance_st;
	TextInput::TimeDataD elevation_st;

	/// Files names for additional environmental input files
	xtring file_disturbance;
	xtring file_disturbance_st;
	xtring file_elevation_st;

	/// Gets disturbance intervals (in years)
	void getdisturbance(Gridcell& gridcell);

	/// Gets elevation
	void getelevation(Gridcell& gridcell);
};

/// Class that deals with all land cover input from text files
class LandcoverInput {

public:

	/// Constructor
	LandcoverInput();

	/// Opens land cover input files.
	void init();

	/// Loads land cover and stand type area fractions from input files
	bool loadlandcover(double lon, double lat);

	/// Gets land cover and stand type fractions for a year.
	/** Updates landcover and stand type variables frac, frac_old and frac_change
	 *  Area fractions are re-scaled if sum is not 1.0
	 */ 
	void getlandcover(Gridcell& gridcell);

	/// Gets crop stand type fractions for a year, called from getlandcover() 
	double get_crop_fractions(Gridcell& gridcell, int year, TimeDataD& CFTdata, double sum_tot);

	/// Gets land cover or stand type transitions for a year
	bool get_land_transitions(Gridcell& gridcell);

	/// Gets land cover transitions for a year
	/** Updates landcover frac_transfer array
	 *  Transition values are checked against net lcc fractions and 
	 *  rescaled if necessary. 
	 */ 
	bool get_lc_transfer(Gridcell& gridcell);

	/// Gets first historic year of net land cover fraction input data
	int getfirsthistyear();

private:

	// Objects handling land cover fraction data input
	TextInput::TimeDataD LUdata;
	TextInput::TimeDataD Peatdata;
	TextInput::TimeDataD grossLUC;
	TextInput::TimeDataD st_data[NLANDCOVERTYPES];

	/// Files names for land cover fraction input files
	xtring file_lu, file_grossLUC, file_peat;
	xtring file_lu_st[NLANDCOVERTYPES];

	/// Whether pfts not in crop fraction input file are removed from pftlist (0,1)
	bool minimizecftlist;

	/// Number of years to increase cropland fraction linearly from 0 to first year's value
	int nyears_cropland_ramp;

	/// Whether to use stand types with suitable rainfed crops (based on crop pft tb and gridcell latitude) when using fixed crop fractions
	bool frac_fixed_default_crops;

	/// Input precision set in instruction file
	int input_precision_force;
};

/// Class that deals with all crop management input from text files
class ManagementInput {

public:

	/// Constructor
	ManagementInput();
	/// Deconstructor
	~ManagementInput();
	/// Opens management data files
	void init();
	/// Loads fertilisation, sowing and harvest dates from input files
	bool loadmanagement(double lon, double lat);
	/// Gets management data for a year
	void getmanagement(Gridcell& gridcell, LandcoverInput& landcover_input);

private:

	/// Input objects for each management text input file
	TextInput::TimeDataD sdates;
	TextInput::TimeDataD hdates;
	TextInput::TimeDataD Nfert;
	TextInput::TimeDataD Nfert_st;
	TextInput::TimeDataD NfertMan;
	TextInput::TimeDataD woodharv_frac;
	TextInput::TimeDataD woodharv_cmass;
	TextInput::TimeDataD cutinterval_st;
	TextInput::TimeDataD firstmanageyear_st;

	TextInput::TimeDataD* targetfrac_pft_mt;

	/// Files names for management input file
	xtring file_sdates, file_hdates, file_Nfert, file_Nfert_st, file_NfertMan, file_woodharv_frac, file_woodharv_cmass,
		file_cutinterval_st, file_firstmanageyear_st;

	/// Gets sowing date data for a year
	void getsowingdates(Gridcell& gridcell);
	/// Gets harvest date data for a year
	void getharvestdates(Gridcell& gridcell);
	/// Gets nitrogen fertilisation data for a year
	void getNfert(Gridcell& gridcell);
	/// Gets wood harvest data for a year
	void getwoodharvest(Gridcell& gridcell, LandcoverInput& landcover_input);
	/// Gets cutinterval
	void getcutinterval(Gridcell& gridcell);
	/// Gets fistmanageyear
	void getfirstmanageyear(Gridcell& gridcell);
	/// Gets pft selection target fractions (per st)
	void gettargetcutting(Gridcell& gridcell);
};

#endif // LPJ_GUESS_EXTERNALINPUT_H