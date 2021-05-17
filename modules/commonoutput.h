///////////////////////////////////////////////////////////////////////////////////////
/// \file commonoutput.h
/// \brief Output module for the most commonly needed output files
///
/// \author Joe Siltberg
/// $Date: 2019-10-28 18:48:52 +0100 (Mon, 28 Oct 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_COMMON_OUTPUT_H
#define LPJ_GUESS_COMMON_OUTPUT_H

#include "outputmodule.h"
#include "outputchannel.h"
#include "gutil.h"

namespace GuessOutput {

/// Output module for the most commonly needed output files
class CommonOutput : public OutputModule {
public:

	CommonOutput();

	~CommonOutput();

	// implemented functions inherited from OutputModule
	// (see documentation in OutputModule)

	void init();

	void outannual(Gridcell& gridcell);

	void outdaily(Gridcell& gridcell);

	void openlocalfiles(Gridcell& gridcell) {};

	void closelocalfiles(Gridcell& gridcell) {};

private:

	/// Defines all output tables
	void define_output_tables();

	// Output file names ...
	xtring file_cmass,file_anpp,file_agpp,file_fpc,file_aaet,file_dens,file_lai,file_cflux,file_doc,file_cpool,file_clitter,file_runoff;
	xtring file_mnpp,file_mlai,file_mgpp,file_mra,file_maet,file_mpet,file_mevap,file_mrunoff,file_mintercep,file_mrh;
	xtring file_mnee,file_mwcont_upper,file_mwcont_lower;
	xtring file_firert,file_speciesheights;
	xtring file_wetland_water_added;

	// bvoc
	xtring file_aiso, file_miso, file_amon, file_mmon, file_amon_mt1, file_amon_mt2, file_mmon_mt1, file_mmon_mt2;

	// nitrogen
	xtring file_nmass, file_cton_leaf, file_nsources, file_npool, file_nlitter, file_nuptake, file_vmaxnlim, file_nflux, file_ngases;
	xtring file_soil_npool, file_soil_nflux;
		
	// BLAZE & SIMFIRE
	xtring file_aburned_area_out, file_mburned_area_out;
	xtring file_simfireanalysis_out;
		
	// Soil temperature at 25cm depth
	xtring file_msoiltempdepth5, file_msoiltempdepth15, file_msoiltempdepth25, file_msoiltempdepth35, file_msoiltempdepth45, file_msoiltempdepth55, file_msoiltempdepth65, file_msoiltempdepth75, file_msoiltempdepth85, file_msoiltempdepth95, file_msoiltempdepth105, file_msoiltempdepth115, file_msoiltempdepth125, file_msoiltempdepth135, file_msoiltempdepth145;
	
	// Methane fluxes
	xtring file_mch4, file_mch4diff, file_mch4plan, file_mch4ebull; 
	
	// Snow, water table depth and active layer depth
	xtring file_msnow, file_mwtp, file_mald;

	// Output tables
	Table out_cmass, out_anpp, out_agpp, out_fpc, out_aaet, out_dens, out_lai, out_cflux, out_doc, out_cpool, out_clitter, out_firert, out_runoff, out_speciesheights;
	Table out_wetland_water_added;

	Table out_mnpp, out_mlai, out_mgpp, out_mra, out_maet, out_mpet, out_mevap, out_mrunoff, out_mintercep;
	Table out_mrh, out_mnee, out_mwcont_upper, out_mwcont_lower;
	
	// bvoc
	Table out_aiso, out_miso, out_amon, out_mmon, out_amon_mt1, out_amon_mt2, out_mmon_mt1, out_mmon_mt2;
	
	Table out_nmass, out_cton_leaf, out_nsources, out_npool, out_nlitter, out_nuptake, out_vmaxnlim, out_nflux, out_ngases;
	Table out_soil_npool, out_soil_nflux;

	// BLAZE && SIMFIRE
	Table out_aburned_area, out_mburned_area;
	Table out_simfireanalysis;

	// Methane, snow, water table and active layer depth
	Table out_mch4, out_mch4diff, out_mch4plan, out_mch4ebull, out_msnow, out_mwtp, out_mald;

	// Soil temperatures
	Table out_msoiltempdepth5, out_msoiltempdepth15, out_msoiltempdepth25, out_msoiltempdepth35, out_msoiltempdepth45, out_msoiltempdepth55, out_msoiltempdepth65, out_msoiltempdepth75, out_msoiltempdepth85, out_msoiltempdepth95, out_msoiltempdepth105, out_msoiltempdepth115, out_msoiltempdepth125, out_msoiltempdepth135, out_msoiltempdepth145;
};

}

#endif // LPJ_GUESS_COMMON_OUTPUT_H
