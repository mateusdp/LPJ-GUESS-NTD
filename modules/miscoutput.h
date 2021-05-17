///////////////////////////////////////////////////////////////////////////////////////
/// \file cropoutput.h
/// \brief Output module for the land use and management information
///
/// \author Joe Siltberg
/// $Date: 2020-09-16 19:22:00 +0200 (Wed, 16 Sep 2020) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_MISC_OUTPUT_H
#define LPJ_GUESS_MISC_OUTPUT_H

#include "outputmodule.h"
#include "outputchannel.h"
#include "gutil.h"

namespace GuessOutput {

// Definitions for separate output files per NATURAL and FOREST stand (when
// instruction file parameter printseparatestands == true, printseparatestands
// is set to false in LandcoverInput::init() when input land cover
// fraction data file has data for > 50 gridcells)

/// Output module for the most commonly needed output files
class MiscOutput : public OutputModule {
public:

	MiscOutput();

	~MiscOutput();

	// implemented functions inherited from OutputModule
	// (see documentation in OutputModule)

	void init();

	void outannual(Gridcell& gridcell);

	void outdaily(Gridcell& gridcell);

private:

	/// Upper limit for files in multiple stand printout
	static const int MAXNUMBER_STANDS = 1000;

	/// Printout of first stand from first historic year
	static const bool PRINTFIRSTSTANDFROM1901 = true;

	/// Defines all output tables
	void openlocalfiles(Gridcell& gridcell);

	void define_output_tables();

	void closelocalfiles(Gridcell& gridcell);

	// Output file names ...
	xtring file_yield, file_yield1, file_yield2, file_sdate1, file_sdate2,
		   file_hdate1, file_hdate2, file_lgp, file_phu, file_fphu, file_fhi,
		   file_irrigation, file_seasonality, file_cflux_cropland,
		   file_cflux_pasture, file_cflux_natural, file_cflux_forest,
		   file_cpool_cropland, file_cpool_pasture, file_cpool_natural,
		   file_cpool_forest, file_nflux_cropland, file_nflux_pasture,
		   file_nflux_natural, file_nflux_forest, file_npool_cropland,
		   file_npool_pasture, file_npool_natural, file_npool_forest,
		   file_anpp_cropland, file_anpp_pasture, file_anpp_natural,
		   file_anpp_forest, file_cmass_cropland, file_cmass_pasture,
		   file_cmass_natural, file_cmass_forest, file_dens_natural,
		   file_dens_forest, file_soil_nflux_cropland, file_soil_nflux_pasture,
		   file_soil_nflux_natural, file_soil_nflux_forest,
		   file_cmass_peatland, file_cflux_peatland,
		   file_cpool_peatland, file_nflux_peatland, file_npool_peatland,
		   file_anpp_peatland;

	// daily
	xtring file_daily_lai, file_daily_npp, file_daily_nmass, file_daily_cmass,
		   file_daily_cton, file_daily_ndemand, file_daily_cmass_leaf,
		   file_daily_nmass_leaf, file_daily_cmass_root, file_daily_nmass_root,
		   file_daily_cmass_stem, file_daily_nmass_stem,
		   file_daily_cmass_storage, file_daily_nmass_storage,
		   file_daily_n_input_soil, file_daily_avail_nmass_soil,
		   file_daily_upper_wcont, file_daily_lower_wcont,
		   file_daily_irrigation, file_daily_climate,
		   file_daily_cmass_dead_leaf, file_daily_nmass_dead_leaf, 
		   file_daily_fphu, file_daily_nminleach,
		   file_daily_norgleach, file_daily_nuptake, file_daily_ds,
		   file_daily_stem, file_daily_leaf, file_daily_root,
		   file_daily_storage;

	// Output tables
	Table out_yield, out_yield1, out_yield2, out_sdate1, out_sdate2,
		  out_hdate1, out_hdate2, out_lgp, out_phu, out_fhi, out_fphu,
		  out_irrigation, out_seasonality, out_cflux_cropland,
		  out_cflux_pasture, out_cflux_natural, out_cflux_forest,
		  out_cpool_cropland, out_cpool_pasture, out_cpool_natural,
		  out_cpool_forest, out_nflux_cropland, out_nflux_pasture,
		  out_nflux_natural, out_nflux_forest, out_npool_cropland,
		  out_npool_pasture, out_npool_natural, out_npool_forest,
		  out_anpp_cropland, out_anpp_pasture, out_anpp_natural,
		  out_anpp_forest, out_cmass_cropland, out_cmass_pasture,
		  out_cmass_natural, out_cmass_forest, out_dens_natural,
		  out_dens_forest, out_soil_nflux_cropland, out_soil_nflux_pasture,
		  out_soil_nflux_natural, out_soil_nflux_forest,
		  out_cflux_peatland, out_cpool_peatland,
		  out_nflux_peatland, out_npool_peatland, out_cmass_peatland,
		  out_anpp_peatland;

	Table* out_anpp_stand[MAXNUMBER_STANDS];
	Table* out_cmass_stand[MAXNUMBER_STANDS];

	//daily
	Table out_daily_lai, out_daily_npp, out_daily_cton, out_daily_nmass,
		  out_daily_cmass, out_daily_ndemand, out_daily_cmass_leaf,
		  out_daily_nmass_leaf, out_daily_cmass_root, out_daily_nmass_root,
		  out_daily_cmass_stem, out_daily_nmass_stem, out_daily_cmass_storage,
		  out_daily_nmass_storage, out_daily_n_input_soil,
		  out_daily_cmass_dead_leaf, out_daily_nmass_dead_leaf, out_daily_fphu,
		  out_daily_avail_nmass_soil, out_daily_upper_wcont,
		  out_daily_lower_wcont, out_daily_irrigation, out_daily_climate,
		  out_daily_nminleach, out_daily_norgleach, out_daily_nuptake, out_daily_ds, 
		  out_daily_stem, out_daily_leaf, out_daily_root, out_daily_storage;
};

}

#endif // LPJ_GUESS_MISC_OUTPUT_H
