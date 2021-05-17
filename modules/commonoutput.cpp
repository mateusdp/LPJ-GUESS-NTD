///////////////////////////////////////////////////////////////////////////////////////
/// \file outputmodule.cpp
/// \brief Implementation of the common output module
///
/// \author Joe Siltberg
/// $Date: 2019-10-28 18:48:52 +0100 (Mon, 28 Oct 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "commonoutput.h"
#include "parameters.h"
#include "guess.h"

// Years between output time points in plot function output (Windows shell only)
const int PLOT_INTERVAL = 1;

// Years between updates of soil water and stand structure plot (Windows shell only)
const int PLOT_UPDATE_INTERVAL = 20;

// Years between updates of 3D vegetation view (Windows shell only)
const int VEG3D_UPDATE_INTERVAL = 5;

namespace GuessOutput {

REGISTER_OUTPUT_MODULE("common", CommonOutput)

CommonOutput::CommonOutput() {
	// Annual output variables
	declare_parameter("file_cmass", &file_cmass, 300, "C biomass output file");
	declare_parameter("file_anpp", &file_anpp, 300, "Annual NPP output file");
	declare_parameter("file_agpp", &file_agpp, 300, "Annual GPP output file");
	declare_parameter("file_fpc", &file_fpc, 300, "FPC output file");
	declare_parameter("file_aaet", &file_aaet, 300, "Annual AET output file");
	declare_parameter("file_lai", &file_lai, 300, "LAI output file");
	declare_parameter("file_cflux", &file_cflux, 300, "C fluxes output file");
	declare_parameter("file_doc", &file_doc, 300, "DOC output file");
	declare_parameter("file_dens", &file_dens, 300, "Tree density output file");
	declare_parameter("file_cpool", &file_cpool, 300, "Soil C output file");
	declare_parameter("file_clitter", &file_clitter, 300, "Litter C output file");
	declare_parameter("file_runoff", &file_runoff, 300, "Runoff output file");
	declare_parameter("file_wetland_water_added", &file_wetland_water_added, 300, "Wetland water added output file");

	declare_parameter("file_firert", &file_firert, 300, "Fire retrun time output file");

	declare_parameter("file_nmass", &file_nmass, 300, "N biomass output file");
	declare_parameter("file_cton_leaf", &file_cton_leaf, 300, "Mean leaf C:N output file");
	declare_parameter("file_nsources", &file_nsources, 300, "Annual nitrogen sources output file");
	declare_parameter("file_npool", &file_npool, 300, "Soil nitrogen output file");
	declare_parameter("file_nlitter", &file_nlitter, 300, "Litter nitrogen output file");
	declare_parameter("file_nuptake", &file_nuptake, 300, "Annual nitrogen uptake output file");
	declare_parameter("file_vmaxnlim", &file_vmaxnlim, 300, "Annual nitrogen limitation on vm output file");
	declare_parameter("file_nflux", &file_nflux, 300, "Annual nitrogen fluxes output file");
	declare_parameter("file_ngases", &file_ngases, 300, "Annual nitrogen gases output file");
	declare_parameter("file_soil_npool", &file_soil_npool, 300, "Annual soil N pools output file");
	declare_parameter("file_soil_nflux", &file_soil_nflux, 300, "Annual soil N fluxes output file");

	declare_parameter("file_speciesheights", &file_speciesheights, 300, "Mean species heights");

	// Monthly output variables
	declare_parameter("file_mnpp", &file_mnpp, 300, "Monthly NPP output file");
	declare_parameter("file_mlai", &file_mlai, 300, "Monthly LAI output file");
	declare_parameter("file_mgpp", &file_mgpp, 300, "Monthly GPP-LeafResp output file");
	declare_parameter("file_mra", &file_mra, 300, "Monthly autotrophic respiration output file");
	declare_parameter("file_maet", &file_maet, 300, "Monthly AET output file");
	declare_parameter("file_mpet", &file_mpet, 300, "Monthly PET output file");
	declare_parameter("file_mevap", &file_mevap, 300, "Monthly Evap output file");
	declare_parameter("file_mrunoff", &file_mrunoff, 300, "Monthly runoff output file");
	declare_parameter("file_mintercep", &file_mintercep, 300, "Monthly intercep output file");
	declare_parameter("file_mrh", &file_mrh, 300, "Monthly heterotrophic respiration output file");
	declare_parameter("file_mnee", &file_mnee, 300, "Monthly NEE output file");
	declare_parameter("file_mwcont_upper", &file_mwcont_upper, 300, "Monthly wcont_upper output file");
	declare_parameter("file_mwcont_lower", &file_mwcont_lower, 300, "Monthly wcont_lower output file");
	// bvoc
	declare_parameter("file_aiso", &file_aiso, 300, "annual isoprene flux output file");
	declare_parameter("file_miso", &file_miso, 300, "monthly isoprene flux output file");
	declare_parameter("file_amon", &file_amon, 300, "annual monoterpene flux output file");
	declare_parameter("file_mmon", &file_mmon, 300, "monthly monoterpene flux output file");
	declare_parameter("file_amon_mt1", &file_amon_mt1, 300, "annual endocyclic monoterpene flux output file");	
	declare_parameter("file_amon_mt2", &file_amon_mt2, 300, "annual other monoterpene flux output file");
	declare_parameter("file_mmon_mt1", &file_mmon_mt1, 300, "monthly endocyclic monoterpene flux output file");	
	declare_parameter("file_mmon_mt2", &file_mmon_mt2, 300, "monthly other monoterpene flux output file");

	if ( firemodel == BLAZE ) {
		declare_parameter("file_aburned_area_out", &file_aburned_area_out, 300, "BLAZE burned area output file");
		declare_parameter("file_mburned_area_out", &file_mburned_area_out, 300, "BLAZE monthly burned area output file");
		declare_parameter("file_simfireanalysis_out", &file_simfireanalysis_out, 300, "SIMFIRE analytics output");
	}

	declare_parameter("file_msoiltempdepth5", &file_msoiltempdepth5, 300, "Soil temperature output file (5cm depth)");
	declare_parameter("file_msoiltempdepth15", &file_msoiltempdepth15, 300, "Soil temperature output file (15cm depth)");
	declare_parameter("file_msoiltempdepth25", &file_msoiltempdepth25, 300, "Soil temperature output file (25cm depth)");
	declare_parameter("file_msoiltempdepth35", &file_msoiltempdepth35, 300, "Soil temperature output file (35cm depth)");
	declare_parameter("file_msoiltempdepth45", &file_msoiltempdepth45, 300, "Soil temperature output file (45cm depth)");
	declare_parameter("file_msoiltempdepth55", &file_msoiltempdepth55, 300, "Soil temperature output file (55cm depth)");
	declare_parameter("file_msoiltempdepth65", &file_msoiltempdepth65, 300, "Soil temperature output file (65cm depth)");
	declare_parameter("file_msoiltempdepth75", &file_msoiltempdepth75, 300, "Soil temperature output file (75cm depth)");
	declare_parameter("file_msoiltempdepth85", &file_msoiltempdepth85, 300, "Soil temperature output file (85cm depth)");
	declare_parameter("file_msoiltempdepth95", &file_msoiltempdepth95, 300, "Soil temperature output file (95cm depth)");
	declare_parameter("file_msoiltempdepth105", &file_msoiltempdepth105, 300, "Soil temperature output file (105cm depth)");
	declare_parameter("file_msoiltempdepth115", &file_msoiltempdepth115, 300, "Soil temperature output file (115cm depth)");
	declare_parameter("file_msoiltempdepth125", &file_msoiltempdepth125, 300, "Soil temperature output file (125cm depth)");
	declare_parameter("file_msoiltempdepth135", &file_msoiltempdepth135, 300, "Soil temperature output file (135cm depth)");
	declare_parameter("file_msoiltempdepth145", &file_msoiltempdepth145, 300, "Soil temperature output file (145cm depth)");

	declare_parameter("file_mch4", &file_mch4, 300, "Monthly CH4 emissions, total");
	declare_parameter("file_mch4diff", &file_mch4diff, 300, "Monthly CH4 emissions, diffusion");
	declare_parameter("file_mch4plan", &file_mch4plan, 300, "Monthly CH4 emissions, plant-mediated");
	declare_parameter("file_mch4ebull", &file_mch4ebull, 300, "Monthly CH4 emissions, ebullition");
	declare_parameter("file_msnow", &file_msnow, 300, "Monthly snow depth");
	declare_parameter("file_mwtp", &file_mwtp, 300, "Monthly water table depth");
	declare_parameter("file_mald", &file_mald, 300, "Monthly active layer depth");
}


CommonOutput::~CommonOutput() {
}

/// Define all output tables and their formats
void CommonOutput::init() {

	define_output_tables();
}

/** This function specifies all columns in all output tables, their names,
 *  column widths and precision.
 *
 *  For each table a TableDescriptor object is created which is then sent to
 *  the output channel to create the table.
 */
void CommonOutput::define_output_tables() {

	//Extra number of decimals when output for benchmarks
#ifdef RUN_BENCHMARKS	
	const int bm_extra_prec = 2;
#else
	const int bm_extra_prec = 0;
#endif
	
	// create a vector with the pft names
	std::vector<std::string> pfts;

	// create a vector with the crop pft names
	std::vector<std::string> crop_pfts;

	pftlist.firstobj();
	while (pftlist.isobj) {
		 Pft& pft=pftlist.getobj();

		 pfts.push_back((char*)pft.name);

		 if(pft.landcover==CROPLAND)
			 crop_pfts.push_back((char*)pft.name);

		 pftlist.nextobj();
	}

	// create a vector with the landcover column titles
	std::vector<std::string> landcovers;

	if (run_landcover) {
		 const char* landcover_string[]={"Urban_sum", "Crop_sum", "Pasture_sum",
			 "Forest_sum", "Natural_sum", "Peatland_sum", "Barren_sum"};
		 for (int i=0; i<NLANDCOVERTYPES; i++) {
			  if(run[i]) {
					landcovers.push_back(landcover_string[i]);
			  }
		 }
	}

	// Create the month columns
	ColumnDescriptors month_columns;
	ColumnDescriptors month_columns_wide;
	xtring months[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
	for (int i = 0; i < 12; i++) {
		month_columns      += ColumnDescriptor(months[i], 8,  3);
		month_columns_wide += ColumnDescriptor(months[i], 10, 3);
	}

	// Create the columns for each output file

	// CMASS
	ColumnDescriptors cmass_columns;
	cmass_columns += ColumnDescriptors(pfts,               8, 3);
	cmass_columns += ColumnDescriptor("Total",             8, 3);
	cmass_columns += ColumnDescriptors(landcovers,        13, 3);

	// ALD
	ColumnDescriptors mald_columns = month_columns;
	mald_columns += ColumnDescriptor("MAXALD", 8, 3);

	// ANPP
	ColumnDescriptors anpp_columns = cmass_columns;

	// AGPP
	ColumnDescriptors agpp_columns = cmass_columns;

	// FPC
	ColumnDescriptors fpc_columns = cmass_columns;

	// AET
	ColumnDescriptors aaet_columns;
	aaet_columns += ColumnDescriptors(pfts,                8, 2);
	aaet_columns += ColumnDescriptor("Total",              8, 2);
	aaet_columns += ColumnDescriptors(landcovers,         13, 2);

	// DENS
	ColumnDescriptors dens_columns;
	dens_columns += ColumnDescriptors(pfts,                8, 4);
	dens_columns += ColumnDescriptor("Total",              8, 4);
	dens_columns += ColumnDescriptors(landcovers,         13, 4);

	// LAI
	ColumnDescriptors lai_columns = dens_columns;

	// CFLUX
	ColumnDescriptors cflux_columns;
	cflux_columns += ColumnDescriptor("Veg",               8, 3);
	cflux_columns += ColumnDescriptor("Repr",              8, 3);
	cflux_columns += ColumnDescriptor("Soil",              8, 3);
	cflux_columns += ColumnDescriptor("Fire",             10, 5);
	cflux_columns += ColumnDescriptor("Est",               8, 3);
	if (run_landcover) {
		 cflux_columns += ColumnDescriptor("Seed",         8, 3);
		 cflux_columns += ColumnDescriptor("Harvest",      9, 3);
		 cflux_columns += ColumnDescriptor("LU_ch",        9, 3);
		 cflux_columns += ColumnDescriptor("Slow_h",       9, 3);
	}
	cflux_columns += ColumnDescriptor("NEE",              10 + bm_extra_prec, 5 + bm_extra_prec);

	ColumnDescriptors doc_columns;
	doc_columns += ColumnDescriptor("Total",              10, 3);
	doc_columns += ColumnDescriptors(landcovers,          13, 3);

	// CPOOL
	ColumnDescriptors cpool_columns;
	cpool_columns += ColumnDescriptor("VegC",              8, 3);

	if (!ifcentury) {
		cpool_columns += ColumnDescriptor("LittC",         8, 3);
		cpool_columns += ColumnDescriptor("SoilfC",        8, 3);
		cpool_columns += ColumnDescriptor("SoilsC",        8, 3);
	}
	else {
		cpool_columns += ColumnDescriptor("LitterC",       8, 3);
		cpool_columns += ColumnDescriptor("SoilC",         8, 3);
	}
	if (run_landcover && ifslowharvestpool) {
		 cpool_columns += ColumnDescriptor("HarvSlowC",   10, 3);
	}
	cpool_columns += ColumnDescriptor("Total",            10 + bm_extra_prec, 3 + bm_extra_prec);

	// CLITTER
	ColumnDescriptors clitter_columns = cmass_columns;

	// FIRERT
	ColumnDescriptors firert_columns;
	firert_columns += ColumnDescriptor("FireRT",			8, 1);
	firert_columns += ColumnDescriptor("BurntAr",			8, 5);

	// BLAZE burnt area 
	ColumnDescriptors blaze_columns;
	blaze_columns += ColumnDescriptor("BurntAr",			9, 5);

	// SIMFIRE Analysis 
	ColumnDescriptors simfireanalysis_columns;
	simfireanalysis_columns += ColumnDescriptor("Biome",	6, 0);
	simfireanalysis_columns += ColumnDescriptor("MxNest",	7, 0);
	simfireanalysis_columns += ColumnDescriptor("PopDens",	10, 3);
	simfireanalysis_columns += ColumnDescriptor("Region",	7, 0);
	
	// RUNOFF
	ColumnDescriptors runoff_columns;
	runoff_columns += ColumnDescriptor("Surf",             8, 1);
	runoff_columns += ColumnDescriptor("Drain",            8, 1);
	runoff_columns += ColumnDescriptor("Base",             8, 1);
	runoff_columns += ColumnDescriptor("Total",            9, 1);

	// WETLAND WATER ADDED
	ColumnDescriptors wetland_water_added_columns;
	wetland_water_added_columns += ColumnDescriptor("H2OAdded", 10, 1);

	// SPECIESHEIGHTS
	ColumnDescriptors speciesheights_columns;
	speciesheights_columns += ColumnDescriptors(pfts,      8, 2);

	// AISO
	ColumnDescriptors aiso_columns;
	aiso_columns += ColumnDescriptors(pfts,               10, 3);
	aiso_columns += ColumnDescriptor("Total",             10, 3);
	aiso_columns += ColumnDescriptors(landcovers,         13, 3);

	// MONOTERPENES
	ColumnDescriptors amon_columns = aiso_columns;

	// CTON
	ColumnDescriptors cton_columns;
	cton_columns += ColumnDescriptors(pfts,                8, 1);
	cton_columns += ColumnDescriptor("Total",              8, 1);
	cton_columns += ColumnDescriptors(landcovers,         13, 1);

	// NSOURCES
	ColumnDescriptors nsources_columns;
	nsources_columns += ColumnDescriptor("NH4dep",         8, 2);
	nsources_columns += ColumnDescriptor("NO3dep",         8, 2);
	nsources_columns += ColumnDescriptor("fix",            8, 2);
	nsources_columns += ColumnDescriptor("fert",           8, 2);
	nsources_columns += ColumnDescriptor("input",          8, 2);
	nsources_columns += ColumnDescriptor("min",            8, 2);
	nsources_columns += ColumnDescriptor("imm",            8, 2);
	nsources_columns += ColumnDescriptor("netmin",         8, 2);
	nsources_columns += ColumnDescriptor("Total",          8, 2);

	// NPOOL
	ColumnDescriptors npool_columns;
	npool_columns += ColumnDescriptor("VegN",              9, 4);
	npool_columns += ColumnDescriptor("LitterN",           9, 4);
	npool_columns += ColumnDescriptor("SoilN",             9, 4);

	if (run_landcover && ifslowharvestpool) {
		npool_columns += ColumnDescriptor("HarvSlowN",    10, 4);
	}

	npool_columns += ColumnDescriptor("Total",            10 + bm_extra_prec, 4 + bm_extra_prec);

	// NMASS
	ColumnDescriptors nmass_columns;
	nmass_columns += ColumnDescriptors(pfts,               8, 2);
	nmass_columns += ColumnDescriptor("Total",             8, 2);
	nmass_columns += ColumnDescriptors(landcovers,        13, 2);

	// NUPTAKE
	ColumnDescriptors nuptake_columns = nmass_columns;

	// NLITTER
	ColumnDescriptors nlitter_columns = nmass_columns;

	// VMAXNLIM
	ColumnDescriptors vmaxnlim_columns;
	vmaxnlim_columns += ColumnDescriptors(pfts,            8, 2);
	vmaxnlim_columns += ColumnDescriptor("Total",          8, 2);
	vmaxnlim_columns += ColumnDescriptors(landcovers,     13, 2);

	// NFLUX
	ColumnDescriptors nflux_columns;
	nflux_columns += ColumnDescriptor("NH4dep",            8, 2);
	nflux_columns += ColumnDescriptor("NO3dep",            8, 2);
	nflux_columns += ColumnDescriptor("fix",               8, 2);
	nflux_columns += ColumnDescriptor("fert",              8, 2);
	nflux_columns += ColumnDescriptor("flux",              8, 2);
	nflux_columns += ColumnDescriptor("leach",             8, 2);
	if (run_landcover) {
		nflux_columns += ColumnDescriptor("seed",		   8, 2);
		nflux_columns += ColumnDescriptor("harvest",       8, 2);
		nflux_columns += ColumnDescriptor("LU_ch",         8, 3);
		nflux_columns += ColumnDescriptor("Slow_h",        8, 3);
	}
	nflux_columns += ColumnDescriptor("NEE",               8 + bm_extra_prec, 2 + bm_extra_prec);

	// NGASES
	ColumnDescriptors ngases_columns;
	ngases_columns += ColumnDescriptor("NH3_fire",         9, 4);
	ngases_columns += ColumnDescriptor("NH3_soil",         9, 4);
	ngases_columns += ColumnDescriptor("NOx_fire",          9, 4);
	ngases_columns += ColumnDescriptor("NOx_soil",          9, 4);
	ngases_columns += ColumnDescriptor("N2O_fire",         9, 4);
	ngases_columns += ColumnDescriptor("N2O_soil",         9, 4);
	ngases_columns += ColumnDescriptor("N2_fire",          9, 4);
	ngases_columns += ColumnDescriptor("N2_soil",          9, 4);
	ngases_columns += ColumnDescriptor("Total",            9, 4);

	// SOIL N TRANSFORMATION - pools
	ColumnDescriptors soil_npool_columns;
	soil_npool_columns += ColumnDescriptor("NH4", 11, 4);
	soil_npool_columns += ColumnDescriptor("NO3", 11, 4);
	soil_npool_columns += ColumnDescriptor("NO2", 11, 4);
	soil_npool_columns += ColumnDescriptor("NO",  11, 4);
	soil_npool_columns += ColumnDescriptor("N2O", 11, 4);
	soil_npool_columns += ColumnDescriptor("N2",  11, 4);

	
	// SOIL N TRANSFORMATION - fluxes
	ColumnDescriptors soil_nflux_columns;
	soil_nflux_columns += ColumnDescriptor("NH3",  12, 6);
	soil_nflux_columns += ColumnDescriptor("NO",   12, 6);
	soil_nflux_columns += ColumnDescriptor("N2O",  12, 6);
	soil_nflux_columns += ColumnDescriptor("N2",   12, 6);

	// *** ANNUAL OUTPUT VARIABLES ***

	create_output_table(out_cmass,          file_cmass,          cmass_columns);
	create_output_table(out_anpp,           file_anpp,           anpp_columns);
	create_output_table(out_agpp,           file_agpp,           agpp_columns);
	create_output_table(out_fpc,            file_fpc,            fpc_columns);
	create_output_table(out_aaet,           file_aaet,           aaet_columns);
	create_output_table(out_dens,           file_dens,           dens_columns);
	create_output_table(out_lai,            file_lai,            lai_columns);
	create_output_table(out_cflux,          file_cflux,          cflux_columns);
	create_output_table(out_doc,	        file_doc,			 doc_columns);
	create_output_table(out_cpool,          file_cpool,          cpool_columns);
	create_output_table(out_clitter,        file_clitter,        clitter_columns);

	if ( firemodel == BLAZE ) {
		create_output_table(out_aburned_area,		file_aburned_area_out,		blaze_columns);
		create_output_table(out_simfireanalysis,	file_simfireanalysis_out,	simfireanalysis_columns);
	} else if ( firemodel == GLOBFIRM ) {
		create_output_table(out_firert,			file_firert,			firert_columns);
	}

	create_output_table(out_runoff,			file_runoff,         runoff_columns);
	create_output_table(out_wetland_water_added, file_wetland_water_added, wetland_water_added_columns);
	create_output_table(out_speciesheights, file_speciesheights, speciesheights_columns);
	create_output_table(out_aiso,           file_aiso,           aiso_columns);
	create_output_table(out_amon,           file_amon,           amon_columns);
	create_output_table(out_amon_mt1,       file_amon_mt1,       amon_columns);
	create_output_table(out_amon_mt2,       file_amon_mt2,       amon_columns);

	create_output_table(out_nmass,          file_nmass,          nmass_columns);
	create_output_table(out_cton_leaf,      file_cton_leaf,      cton_columns);
	create_output_table(out_nsources,       file_nsources,       nsources_columns);
	create_output_table(out_npool,          file_npool,          npool_columns);
	create_output_table(out_nlitter,        file_nlitter,        nlitter_columns);
	create_output_table(out_nuptake,        file_nuptake,        nuptake_columns);
	create_output_table(out_vmaxnlim,       file_vmaxnlim,       vmaxnlim_columns);
	create_output_table(out_nflux,          file_nflux,          nflux_columns);
	create_output_table(out_ngases,         file_ngases,         ngases_columns);
	create_output_table(out_soil_npool,		file_soil_npool,     soil_npool_columns);
	create_output_table(out_soil_nflux,		file_soil_nflux,     soil_nflux_columns);

	// *** MONTHLY OUTPUT VARIABLES ***

	create_output_table(out_mnpp,           file_mnpp,           month_columns);
	create_output_table(out_mlai,           file_mlai,           month_columns);
	create_output_table(out_mgpp,           file_mgpp,           month_columns);
	create_output_table(out_mra,            file_mra,            month_columns);
	create_output_table(out_maet,           file_maet,           month_columns);
	create_output_table(out_mpet,           file_mpet,           month_columns);
	create_output_table(out_mevap,          file_mevap,          month_columns);
	create_output_table(out_mrunoff,        file_mrunoff,        month_columns_wide);
	create_output_table(out_mintercep,      file_mintercep,      month_columns);
	create_output_table(out_mrh,            file_mrh,            month_columns);
	create_output_table(out_mnee,           file_mnee,           month_columns);
	create_output_table(out_mwcont_upper,   file_mwcont_upper,   month_columns);
	create_output_table(out_mwcont_lower,   file_mwcont_lower,   month_columns);
	create_output_table(out_miso,           file_miso,           month_columns_wide);
	create_output_table(out_mmon,           file_mmon,           month_columns_wide);
	create_output_table(out_mmon_mt1,       file_mmon_mt1,       month_columns_wide);
	create_output_table(out_mmon_mt2,       file_mmon_mt2,       month_columns_wide);
	create_output_table(out_mburned_area,   file_mburned_area_out, month_columns);
    
	// Methane
	create_output_table(out_mch4,           file_mch4,           month_columns);
	create_output_table(out_mch4diff,       file_mch4diff,       month_columns);
	create_output_table(out_mch4plan,       file_mch4plan,       month_columns);
	create_output_table(out_mch4ebull,      file_mch4ebull,      month_columns);
    
	// Snow
	create_output_table(out_msnow,          file_msnow,          month_columns);
	create_output_table(out_mwtp,           file_mwtp,           month_columns);
	create_output_table(out_mald,           file_mald,           mald_columns);

	// Soil temperatures
	create_output_table(out_msoiltempdepth5, file_msoiltempdepth5, month_columns);
	create_output_table(out_msoiltempdepth15, file_msoiltempdepth15, month_columns);
	create_output_table(out_msoiltempdepth25, file_msoiltempdepth25, month_columns);
	create_output_table(out_msoiltempdepth35, file_msoiltempdepth35, month_columns);
	create_output_table(out_msoiltempdepth45, file_msoiltempdepth45, month_columns);
	create_output_table(out_msoiltempdepth55, file_msoiltempdepth55, month_columns);
	create_output_table(out_msoiltempdepth65, file_msoiltempdepth65, month_columns);
	create_output_table(out_msoiltempdepth75, file_msoiltempdepth75, month_columns);
	create_output_table(out_msoiltempdepth85, file_msoiltempdepth85, month_columns);
	create_output_table(out_msoiltempdepth95, file_msoiltempdepth95, month_columns);
	create_output_table(out_msoiltempdepth105, file_msoiltempdepth105, month_columns);
	create_output_table(out_msoiltempdepth115, file_msoiltempdepth115, month_columns);
	create_output_table(out_msoiltempdepth125, file_msoiltempdepth125, month_columns);
	create_output_table(out_msoiltempdepth135, file_msoiltempdepth135, month_columns);
	create_output_table(out_msoiltempdepth145, file_msoiltempdepth145, month_columns);
}

/// Function for producing data file used to communicate information on stand structure
/** for 3D vegetation plot in Windows shell
 */
void output_vegetation(Gridcell& gridcell, Pftlist& pftlist) {
	
	// File for output of 3D vegetation structure (invoked by Windows shell only)
	plot3d_fileopen(); 
	
	if (plot3d_getfilehandle()) {

		int ival, p, npft_tree, npft_grass, npft_total;
		double grasslai;
		const bool FALSCH = false;
		const double rgb[3] = { -1, -1, -1 };
		char pftname[16];

		// Loop through Stands
		Gridcell::iterator gc_itr = gridcell.begin();
		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			npft_tree = npft_grass = 0;
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				if (pft.lifeform == TREE) npft_tree++;
				else if (pft.lifeform == GRASS) npft_grass++;
				pftlist.nextobj();
			}
			npft_total = npft_tree + npft_grass;
			fwrite(&npft_total, sizeof(int), 1, plot3d_getfilehandle());
			fwrite(&npft_tree, sizeof(int), 1, plot3d_getfilehandle());
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				if (pft.lifeform == TREE) {
					sprintf(pftname, "%s", (char*)(pft.name.left(15)));
					fwrite(pftname, sizeof(char), 16, plot3d_getfilehandle());
					fwrite(&FALSCH, sizeof(bool), 1, plot3d_getfilehandle());
					fwrite(&rgb, sizeof(double), 3, plot3d_getfilehandle());
					//fwrite(&pft.ifconifer, sizeof(bool), 1, plot3d_getfilehandle());
					//fwrite(pft.preferredrgb, sizeof(double), 3, plot3d_getfilehandle());
				}
				pftlist.nextobj();
			}
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				if (pft.lifeform == GRASS) {
					sprintf(pftname, "%s", (char*)(pft.name.left(15)));
					fwrite(pftname, sizeof(char), 16, plot3d_getfilehandle());
					fwrite(&rgb, sizeof(double), 3, plot3d_getfilehandle());
					//fwrite(pft.preferredrgb, sizeof(double), 3, plot3d_getfilehandle());
				}
				pftlist.nextobj();
			}
			int npatch = stand.npatch();
			fwrite(&npatch, sizeof(int), 1, plot3d_getfilehandle());
			fwrite(&patcharea, sizeof(double), 1, plot3d_getfilehandle());
			for (p = 0; p<npatch; p++) {
				Patch& patch = stand[p];
				Vegetation& vegetation = patch.vegetation;
				grasslai = 0.0;
				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv = vegetation.getobj();
					if (indiv.pft.lifeform == TREE && indiv.alive) {
						ival = indiv.pft.id;
						fwrite(&ival, sizeof(int), 1, plot3d_getfilehandle());
						ival = indiv.id;
						fwrite(&ival, sizeof(int), 1, plot3d_getfilehandle());
						fwrite(&indiv.densindiv, sizeof(double), 1, plot3d_getfilehandle());
						fwrite(&indiv.height, sizeof(double), 1, plot3d_getfilehandle());
						fwrite(&indiv.crownarea, sizeof(double), 1, plot3d_getfilehandle());
					}
					else if (indiv.pft.lifeform == GRASS) grasslai += indiv.lai;
					vegetation.nextobj();
				}
				ival = -9999;
				fwrite(&ival, sizeof(int), 1, plot3d_getfilehandle()); // indicates no more cohorts in this patch
				if (grasslai<0) grasslai = 0.0;
				fwrite(&grasslai, sizeof(double), 1, plot3d_getfilehandle());
			}
			++gc_itr;
		} //while (gridcell.isobj) 

		plot3d_fileclose();

	} //if(out)

	plot3d();
}

/// Gets stand age structure to argument densindiv of dimensions [npft,nageclass]
/** First call with havedims=false to allocate memory and return nageclass,
  * then with havedims=true to get data into densindiv
  * Calling function is responsible for deallocating memory by:
  * delete[] densindiv;
*/
void get_stand_age_structure(Gridcell& gridcell,double* densindiv,int& nageclass,bool havedims) {

	int p, c;
	double active_fraction;

	if (!(vegmode == COHORT || vegmode == INDIVIDUAL)) {
		nageclass = 0;
		return;
	}

	if (havedims) {
		for (p = 0; p < npft; p++)
			for (c = 0; c < nageclass; c++)
				densindiv[p*nageclass + c] = 0.0;
	}
	else nageclass = 0;

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();

		// Determine area fraction of stands where this pft is active:
		active_fraction = 0.0;
		Gridcell::iterator gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			if (stand.pft[pft.id].active) {
				active_fraction += stand.get_gridcell_fraction();
			}
			++gc_itr;
		}

		// Loop through Stands
		gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			Standpft& standpft = stand.pft[pft.id];
			if (standpft.active) {

				stand.firstobj();

				// Loop through Patches
				while (stand.isobj) {
					Patch& patch = stand.getobj();
					Patchpft& patchpft = patch.pft[pft.id];
					Vegetation& vegetation = patch.vegetation;

					// Loop through individuals/cohorts

					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv = vegetation.getobj();

						if (indiv.id != -1 && indiv.alive) {

							if (indiv.pft.id == pft.id) {

								// Age structure
								c = (int)(indiv.age / estinterval);

								if (havedims) {

									densindiv[pft.id*nageclass + c] += indiv.densindiv / (double)stand.npatch()* stand.get_gridcell_fraction() / active_fraction;

								}
								else if (c > nageclass) nageclass = c;
							}
						}
						vegetation.nextobj();
					}
					stand.nextobj();
				}
			}
			++gc_itr;
		}
		pftlist.nextobj();
	}
}

/// Local analogue of OutputRows::add_value for restricting output
/** Use to restrict output to specified range of years
  * (or other user-specified limitation)
  *
  * If only yearly output between, say 1961 and 1990 is requred, use:
  *  if (date.get_calendar_year() >= 1961 && date.get_calendar_year() <= 1990)
  *  (assuming the input module has set the first calendar year in the date object)	
  */
void outlimit(OutputRows& out, const Table& table, double d) {

	if (date.year >= nyear_spinup)
		out.add_value(table, d);
}

/// Output of simulation results at the end of each year
/** Output of simulation results at the end of each year, or for specific years in
  * the simulation of each stand or grid cell. 
  * This function does not have to provide any information to the framework.
  *
  * Restrict output to specific years in the local helper function outlimit().
  *
  * Changes in the structure of this function should be mirrored in outannual()
  * of the other output modules, e.g. MiscOutput::outannual().
  */
void CommonOutput::outannual(Gridcell& gridcell) {

	int c, m;
	double flux_veg, flux_repr, flux_soil, flux_fire, flux_est, flux_seed, flux_charvest;
	double c_fast, c_slow, c_harv_slow;

	double surfsoillitterc,surfsoillittern,cwdc,cwdn,centuryc,centuryn,n_harv_slow,availn;
	double flux_nharvest, flux_nseed;
	double flux_NH3_soil,flux_NOx_soil,flux_N2O_soil,flux_N2_soil;
	double flux_NH3_fire,flux_NOx_fire,flux_N2O_fire,flux_N2_fire;
	double flux_ntot;

	double NH4_mass, NH3_mass, NO3_mass;
	double NO2_mass, NO_mass, N2O_mass, N2_mass;
	double NO_mass_inc,N2O_mass_inc,N2_mass_inc,NH3_mass_inc;
	double gross_nitrif, gross_denitrif, net_nitrif, net_denitrif;

	// hold the monthly average across patches
	double mnpp[12];
	double mgpp[12];
	double mlai[12];
	double maet[12];
	double mpet[12];
	double mevap[12];
	double mintercep[12];
	double mrunoff[12];
	double mrh[12];
	double mra[12];
	double mnee[12];
	double mwcont_upper[12];
	double mwcont_lower[12];
	double miso[12];
	double mmon[12];
	double mmon_mt1[12];
	double mmon_mt2[12];

	double msoilt[12][SOILTEMPOUT];
	double mch4[12];
	double mch4_diff[12];
	double mch4_ebull[12];
	double mch4_plant[12];
	double msnowdepth[12];
	double mwtp[12];
	double mald[12];

	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();

	// The OutputRows object manages the next row of output for each
	// output table
	OutputRows out(output_channel, lon, lat, date.get_calendar_year());

	// guess2008 - reset monthly and annual sums across patches each year
	for (m = 0; m < 12; m++) {
		mnpp[m] = mlai[m] = mgpp[m] = mra[m] = maet[m] = mpet[m] = mevap[m] = mintercep[m] = mrunoff[m] = mrh[m] = mnee[m] = mwcont_upper[m] = mwcont_lower[m] = miso[m] = mmon[m] = mmon_mt1[m] = mmon_mt2[m] = 0.0;

		for (int sl = 0; sl < SOILTEMPOUT; sl++) msoilt[m][sl] = 0.0;
		mch4[m] = mch4_diff[m] = mch4_ebull[m] = mch4_plant[m] = msnowdepth[m] = mwtp[m] = mald[m] = 0.0;
	}

	double aaet, apet, aevap, arunoff, aintercep, awetland_water_added;
	aaet = apet = aevap = arunoff = aintercep = awetland_water_added = 0.0;

	double landcover_cmass[NLANDCOVERTYPES]={0.0};
	double landcover_nmass[NLANDCOVERTYPES]={0.0};
	double landcover_cmass_leaf[NLANDCOVERTYPES]={0.0};
	double landcover_nmass_leaf[NLANDCOVERTYPES]={0.0};
	double landcover_cmass_veg[NLANDCOVERTYPES]={0.0};
	double landcover_nmass_veg[NLANDCOVERTYPES]={0.0};
	double landcover_clitter[NLANDCOVERTYPES]={0.0};
	double landcover_nlitter[NLANDCOVERTYPES]={0.0};
	double landcover_anpp[NLANDCOVERTYPES]={0.0};
	double landcover_agpp[NLANDCOVERTYPES]={0.0};
	double landcover_fpc[NLANDCOVERTYPES]={0.0};
	double landcover_aaet[NLANDCOVERTYPES]={0.0};
	double landcover_lai[NLANDCOVERTYPES]={0.0};
	double landcover_densindiv_total[NLANDCOVERTYPES]={0.0};
	double landcover_aiso[NLANDCOVERTYPES]={0.0};
	double landcover_amon[NLANDCOVERTYPES]={0.0};
	double landcover_amon_mt1[NLANDCOVERTYPES]={0.0};
	double landcover_amon_mt2[NLANDCOVERTYPES]={0.0};
	double landcover_nuptake[NLANDCOVERTYPES]={0.0};
	double landcover_vmaxnlim[NLANDCOVERTYPES]={0.0};

	double mean_standpft_cmass=0.0;
	double mean_standpft_nmass=0.0;
	double mean_standpft_cmass_leaf=0.0;
	double mean_standpft_nmass_leaf=0.0;
	double mean_standpft_cmass_veg=0.0;
	double mean_standpft_nmass_veg=0.0;
	double mean_standpft_clitter=0.0;
	double mean_standpft_nlitter=0.0;
	double mean_standpft_anpp=0.0;
	double mean_standpft_agpp=0.0;
	double mean_standpft_fpc=0.0;
	double mean_standpft_aaet=0.0;
	double mean_standpft_lai=0.0;
	double mean_standpft_densindiv_total=0.0;
	double mean_standpft_heightindiv_total=0.0;
	double mean_standpft_aiso=0.0;
	double mean_standpft_amon=0.0;
	double mean_standpft_amon_mt1=0.0;
	double mean_standpft_amon_mt2=0.0;
	double mean_standpft_nuptake=0.0;
	double mean_standpft_vmaxnlim=0.0;

	double cmass_gridcell=0.0;
	double nmass_gridcell= 0.0;
	double cmass_leaf_gridcell=0.0;
	double nmass_leaf_gridcell=0.0;
	double cmass_veg_gridcell=0.0;
	double nmass_veg_gridcell=0.0;
	double clitter_gridcell=0.0;
	double nlitter_gridcell= 0.0;
	double anpp_gridcell=0.0;
	double agpp_gridcell=0.0;
	double fpc_gridcell=0.0;
	double aaet_gridcell=0.0;
	double lai_gridcell=0.0;
	double surfrunoff_gridcell=0.0;
	double drainrunoff_gridcell=0.0;
	double baserunoff_gridcell=0.0;
	double runoff_gridcell=0.0;
	double wetland_water_added_gridcell = 0.0;
	double dens_gridcell=0.0;
	double firert_gridcell=0.0;
	double burned_area_gridcell=0.0;
	double aiso_gridcell=0.0;
	double amon_gridcell=0.0;
	double amon_mt1_gridcell=0.0;
	double amon_mt2_gridcell=0.0;
	double nuptake_gridcell=0.0;
	double vmaxnlim_gridcell=0.0;
	double maxald_gridcell=0.0;

	double aNH4dep_gridcell=0.0;
	double aNO3dep_gridcell=0.0;
	double anfert_gridcell=0.0;
	double anmin_gridcell=0.0;
	double animm_gridcell=0.0;
	double anfix_gridcell=0.0;
	double n_min_leach_gridcell=0.0;
	double n_org_leach_gridcell=0.0;
	double c_org_leach_gridcell=0.0;

	double standpft_cmass=0.0;
	double standpft_nmass=0.0;
	double standpft_cmass_leaf=0.0;
	double standpft_nmass_leaf=0.0;
	double standpft_cmass_veg=0.0;
	double standpft_nmass_veg=0.0;
	double standpft_clitter=0.0;
	double standpft_nlitter=0.0;
	double standpft_anpp=0.0;
	double standpft_agpp=0.0;
	double standpft_fpc=0.0;
	double standpft_aaet=0.0;
	double standpft_lai=0.0;
	double standpft_densindiv_total=0.0;
	double standpft_heightindiv_total = 0.0;
	double standpft_aiso=0.0;
	double standpft_amon=0.0;
	double standpft_amon_mt1=0.0;
	double standpft_amon_mt2=0.0;
	double standpft_nuptake=0.0;
	double standpft_vmaxnlim=0.0;

	// *** Loop through PFTs ***

	pftlist.firstobj();
	while (pftlist.isobj) {

		Pft& pft=pftlist.getobj();
		Gridcellpft& gridcellpft=gridcell.pft[pft.id];

		// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
		mean_standpft_cmass=0.0;
		mean_standpft_nmass=0.0;
		mean_standpft_cmass_leaf=0.0;
		mean_standpft_nmass_leaf=0.0;
		mean_standpft_cmass_veg=0.0;
		mean_standpft_nmass_veg=0.0;
		mean_standpft_clitter=0.0;
		mean_standpft_nlitter=0.0;
		mean_standpft_anpp=0.0;
		mean_standpft_agpp=0.0;
		mean_standpft_fpc=0.0;
		mean_standpft_aaet=0.0;
		mean_standpft_lai=0.0;
		mean_standpft_densindiv_total=0.0;
		mean_standpft_aiso=0.0;
		mean_standpft_amon=0.0;
		mean_standpft_amon_mt1=0.0;
		mean_standpft_amon_mt2=0.0;
		mean_standpft_nuptake=0.0;
		mean_standpft_vmaxnlim=0.0;

		mean_standpft_heightindiv_total = 0.0;

		// Determine area fraction of stands where this pft is active:
		double active_fraction = 0.0;

		Gridcell::iterator gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			if(stand.pft[pft.id].active) {
				active_fraction += stand.get_gridcell_fraction();
			}

			++gc_itr;
		}

		// Loop through Stands
		gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			Standpft& standpft=stand.pft[pft.id];
			if(standpft.active) {
			// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
			standpft_cmass=0.0;
			standpft_nmass=0.0;
			standpft_cmass_leaf=0.0;
			standpft_nmass_leaf=0.0;
			standpft_cmass_veg=0.0;
			standpft_nmass_veg=0.0;
			standpft_clitter=0.0;
			standpft_nlitter=0.0;
			standpft_anpp=0.0;
			standpft_agpp=0.0;
			standpft_fpc=0.0;
			standpft_aaet=0.0;
			standpft_lai=0.0;
			standpft_densindiv_total = 0.0;
			standpft_heightindiv_total = 0.0;
			standpft_aiso=0.0;
			standpft_amon=0.0;
			standpft_amon_mt1=0.0;
			standpft_amon_mt2=0.0;
			standpft_nuptake=0.0;
			standpft_vmaxnlim=0.0;

			stand.firstobj();

			// Loop through Patches
			while (stand.isobj) {
				Patch& patch = stand.getobj();
				Patchpft& patchpft = patch.pft[pft.id];
				Vegetation& vegetation = patch.vegetation;

				standpft_anpp += patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id);
				standpft_agpp += patch.fluxes.get_annual_flux(Fluxes::GPP, pft.id);
				standpft_aiso += patch.fluxes.get_annual_flux(Fluxes::ISO, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_APIN, pft.id); 					
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_LIMO, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_TRIC, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_BPIN, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_MYRC, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_SABI, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_CAMP, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_TBOC, pft.id);
				standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MT_OTHR, pft.id);
				standpft_amon_mt1 += patch.fluxes.get_annual_flux(Fluxes::MT_APIN, pft.id);
				standpft_amon_mt1 += patch.fluxes.get_annual_flux(Fluxes::MT_LIMO, pft.id);
				standpft_amon_mt1 += patch.fluxes.get_annual_flux(Fluxes::MT_TRIC, pft.id);
				standpft_amon_mt2 += patch.fluxes.get_annual_flux(Fluxes::MT_BPIN, pft.id);
				standpft_amon_mt2 += patch.fluxes.get_annual_flux(Fluxes::MT_MYRC, pft.id);
				standpft_amon_mt2 += patch.fluxes.get_annual_flux(Fluxes::MT_SABI, pft.id);
				standpft_amon_mt2 += patch.fluxes.get_annual_flux(Fluxes::MT_CAMP, pft.id);
				standpft_amon_mt2 += patch.fluxes.get_annual_flux(Fluxes::MT_TBOC, pft.id);
				standpft_amon_mt2 += patch.fluxes.get_annual_flux(Fluxes::MT_OTHR, pft.id);
					
				standpft_clitter += patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr;
				standpft_nlitter += patchpft.nmass_litter_leaf + patchpft.nmass_litter_root + patchpft.nmass_litter_sap + patchpft.nmass_litter_heart;

					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv=vegetation.getobj();

						if (indiv.id!=-1 && indiv.alive) {

							if (indiv.pft.id==pft.id) {

								standpft_cmass_leaf += indiv.cmass_leaf;
								standpft_cmass += indiv.ccont();
								standpft_nmass += indiv.ncont();
								standpft_nmass_leaf += indiv.cmass_leaf / indiv.cton_leaf_aavr;
								standpft_nmass_veg += indiv.nmass_veg;
								standpft_fpc += indiv.fpc;
								standpft_aaet += indiv.aaet;
								standpft_lai += indiv.lai;
								if (pft.lifeform==TREE) {	
									standpft_densindiv_total += indiv.densindiv;
									standpft_heightindiv_total += indiv.height * indiv.densindiv;
								}
								standpft_vmaxnlim += indiv.avmaxnlim * indiv.cmass_leaf;
								standpft_nuptake += indiv.anuptake;

								if(pft.landcover == CROPLAND) {
									standpft_cmass_veg += indiv.cmass_leaf + indiv.cmass_root;
									if(indiv.cropindiv) {
										standpft_cmass_veg += indiv.cropindiv->cmass_ho + indiv.cropindiv->cmass_agpool + indiv.cropindiv->cmass_stem;
										standpft_nmass_leaf += indiv.cropindiv->ynmass_leaf + indiv.cropindiv->ynmass_dead_leaf;
										standpft_nmass_veg += indiv.cropindiv->ynmass_leaf + indiv.cropindiv->ynmass_dead_leaf + indiv.cropindiv->ynmass_root + indiv.cropindiv->ynmass_ho + indiv.cropindiv->ynmass_agpool;
									}
								}
								else {

									standpft_cmass_veg += indiv.cmass_veg;

								}
							}

						} // alive?
						vegetation.nextobj();
					}

					stand.nextobj();
				} // end of patch loop

				standpft_cmass/=(double)stand.npatch();
				standpft_nmass/=(double)stand.npatch();
				standpft_cmass_leaf/=(double)stand.npatch();
				standpft_nmass_leaf/=(double)stand.npatch();
				standpft_cmass_veg/=(double)stand.npatch();
				standpft_nmass_veg/=(double)stand.npatch();
				standpft_clitter/=(double)stand.npatch();
				standpft_nlitter/=(double)stand.npatch();
				standpft_anpp/=(double)stand.npatch();
				standpft_agpp/=(double)stand.npatch();
				standpft_fpc/=(double)stand.npatch();
				standpft_aaet/=(double)stand.npatch();
				standpft_lai/=(double)stand.npatch();
				standpft_densindiv_total/=(double)stand.npatch();
				standpft_aiso/=(double)stand.npatch(); // missing above!
				standpft_amon/=(double)stand.npatch(); // missing above!
				standpft_amon_mt1/=(double)stand.npatch();
				standpft_amon_mt2/=(double)stand.npatch();
				standpft_nuptake/=(double)stand.npatch();
				standpft_vmaxnlim/=(double)stand.npatch();
				standpft_heightindiv_total/=(double)stand.npatch();

				if (!negligible(standpft_cmass_leaf))
					standpft_vmaxnlim /= standpft_cmass_leaf;

				//Update landcover totals
				landcover_cmass[stand.landcover]+=standpft_cmass*stand.get_landcover_fraction();
				landcover_nmass[stand.landcover]+=standpft_nmass*stand.get_landcover_fraction();
				landcover_cmass_leaf[stand.landcover]+=standpft_cmass_leaf*stand.get_landcover_fraction();
				landcover_nmass_leaf[stand.landcover]+=standpft_nmass_leaf*stand.get_landcover_fraction();
				landcover_cmass_veg[stand.landcover]+=standpft_cmass_veg*stand.get_landcover_fraction();
				landcover_nmass_veg[stand.landcover]+=standpft_nmass_veg*stand.get_landcover_fraction();
				landcover_clitter[stand.landcover]+=standpft_clitter*stand.get_landcover_fraction();
				landcover_nlitter[stand.landcover]+=standpft_nlitter*stand.get_landcover_fraction();
				landcover_anpp[stand.landcover]+=standpft_anpp*stand.get_landcover_fraction();
				landcover_agpp[stand.landcover]+=standpft_agpp*stand.get_landcover_fraction();
				if(!pft.isintercropgrass) {
					landcover_fpc[stand.landcover]+=standpft_fpc*stand.get_landcover_fraction();
					landcover_lai[stand.landcover]+=standpft_lai*stand.get_landcover_fraction();
				}
				landcover_aaet[stand.landcover]+=standpft_aaet*stand.get_landcover_fraction();
				landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*stand.get_landcover_fraction();
				landcover_aiso[stand.landcover]+=standpft_aiso*stand.get_landcover_fraction();
				landcover_amon[stand.landcover]+=standpft_amon*stand.get_landcover_fraction();
				landcover_amon_mt1[stand.landcover]+=standpft_amon_mt1*stand.get_landcover_fraction();
				landcover_amon_mt2[stand.landcover]+=standpft_amon_mt2*stand.get_landcover_fraction();
				landcover_nuptake[stand.landcover]+=standpft_nuptake*stand.get_landcover_fraction();
				landcover_vmaxnlim[stand.landcover]+=standpft_vmaxnlim*stand.get_landcover_fraction();

				//Update pft means for active stands
				if(active_fraction) {
					mean_standpft_cmass += standpft_cmass * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_nmass += standpft_nmass * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_cmass_leaf += standpft_cmass_leaf * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_nmass_leaf += standpft_nmass_leaf * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_cmass_veg += standpft_cmass_veg * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_nmass_veg += standpft_nmass_veg * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_clitter += standpft_clitter * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_nlitter += standpft_nlitter * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_anpp += standpft_anpp * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_agpp += standpft_agpp * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_fpc += standpft_fpc * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_aaet += standpft_aaet * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_lai += standpft_lai * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_densindiv_total += standpft_densindiv_total * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_heightindiv_total += standpft_heightindiv_total * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_aiso += standpft_aiso * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_amon += standpft_amon * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_amon_mt1 += standpft_amon_mt1 * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_amon_mt2 += standpft_amon_mt2 * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_nuptake += standpft_nuptake * stand.get_gridcell_fraction() / active_fraction;
					mean_standpft_vmaxnlim += standpft_vmaxnlim * stand.get_gridcell_fraction() / active_fraction;
				}

				//Update stand totals
				stand.anpp += standpft_anpp;
				stand.cmass += standpft_cmass;

				// Update gridcell totals
				double fraction_of_gridcell = stand.get_gridcell_fraction();

				cmass_gridcell+=standpft_cmass*fraction_of_gridcell;
				nmass_gridcell+=standpft_nmass*fraction_of_gridcell;
				cmass_leaf_gridcell+=standpft_cmass_leaf*fraction_of_gridcell;
				nmass_leaf_gridcell+=standpft_nmass_leaf*fraction_of_gridcell;
				cmass_veg_gridcell+=standpft_cmass_veg*fraction_of_gridcell;
				nmass_veg_gridcell+=standpft_nmass_veg*fraction_of_gridcell;
				clitter_gridcell+=standpft_clitter*fraction_of_gridcell;
				nlitter_gridcell+=standpft_nlitter*fraction_of_gridcell;
				anpp_gridcell+=standpft_anpp*fraction_of_gridcell;
				agpp_gridcell+=standpft_agpp*fraction_of_gridcell;
				if(!pft.isintercropgrass) {
					fpc_gridcell+=standpft_fpc*fraction_of_gridcell;
					lai_gridcell+=standpft_lai*fraction_of_gridcell;
				}
				aaet_gridcell+=standpft_aaet*fraction_of_gridcell;
				dens_gridcell+=standpft_densindiv_total*fraction_of_gridcell;
				aiso_gridcell+=standpft_aiso*fraction_of_gridcell;
				amon_gridcell+=standpft_amon*fraction_of_gridcell;
				amon_mt1_gridcell+=standpft_amon_mt1*fraction_of_gridcell;
				amon_mt2_gridcell+=standpft_amon_mt2*fraction_of_gridcell;
				nuptake_gridcell+=standpft_nuptake*fraction_of_gridcell;
				vmaxnlim_gridcell+=standpft_vmaxnlim*standpft_cmass_leaf*fraction_of_gridcell;

				// Graphical output every PLOT_INTERVAL years
				// (Windows shell only - "plot" statements have no effect otherwise)
				if (!(date.year%PLOT_INTERVAL)) {
					plot("C mass [kgC/m2]",pft.name,date.year,mean_standpft_cmass);
					plot("NPP [kgC/m2/yr]",pft.name,date.year,mean_standpft_anpp);
					plot("LAI [m2/m2]",pft.name,date.year,mean_standpft_lai);
					if (pft.lifeform == TREE) plot("Dens [indiv/ha]",pft.name,date.year,mean_standpft_densindiv_total*M2_PER_HA);
					if (mean_standpft_cmass_leaf > 0.0 && ifnlim) {
						plot("Vmax N lim",pft.name,date.year,mean_standpft_vmaxnlim);
						plot("leaf C:N [kgC/kg N]",pft.name,date.year,mean_standpft_cmass_leaf/mean_standpft_nmass_leaf);
					}
				}

			}//if(active)
			++gc_itr;
		}//End of loop through stands

		// Print PFT sums to files

		double standpft_mean_cton_leaf = limited_cton(mean_standpft_cmass_leaf, mean_standpft_nmass_leaf);

		outlimit(out,out_cmass,     mean_standpft_cmass);
		outlimit(out,out_anpp,      mean_standpft_anpp);
		outlimit(out,out_agpp,      mean_standpft_agpp);
		outlimit(out,out_fpc,       mean_standpft_fpc);
		outlimit(out,out_aaet,      mean_standpft_aaet);
		outlimit(out,out_clitter,   mean_standpft_clitter);
		outlimit(out,out_dens,      mean_standpft_densindiv_total);
		outlimit(out,out_lai,       mean_standpft_lai);
		outlimit(out,out_aiso,      mean_standpft_aiso);
		outlimit(out,out_amon,      mean_standpft_amon);
		outlimit(out,out_amon_mt1,  mean_standpft_amon_mt1);
		outlimit(out,out_amon_mt2,  mean_standpft_amon_mt2);
		outlimit(out,out_nmass,     (mean_standpft_nmass + mean_standpft_nlitter) * M2_PER_HA);
		outlimit(out,out_cton_leaf, standpft_mean_cton_leaf);
		outlimit(out,out_vmaxnlim,  mean_standpft_vmaxnlim);
		outlimit(out,out_nuptake,   mean_standpft_nuptake * M2_PER_HA);
		outlimit(out,out_nlitter,   mean_standpft_nlitter * M2_PER_HA);

		// print species heights
		double height = 0.0;
		if (mean_standpft_densindiv_total > 0.0)
			height = mean_standpft_heightindiv_total / mean_standpft_densindiv_total;

		outlimit(out,out_speciesheights, height);

		pftlist.nextobj();

	} // *** End of PFT loop ***

	flux_veg = flux_repr = flux_soil = flux_fire = flux_est = flux_seed = flux_charvest = 0.0;

	// guess2008 - carbon pools
	c_fast = c_slow = c_harv_slow = 0.0;

	surfsoillitterc = surfsoillittern = cwdc = cwdn = centuryc = centuryn = n_harv_slow = availn = 0.0;
	aNH4dep_gridcell = aNO3dep_gridcell = anfert_gridcell = anmin_gridcell = animm_gridcell = anfix_gridcell = 0.0;
	n_org_leach_gridcell = n_min_leach_gridcell = c_org_leach_gridcell = 0.0;
	flux_NH3_soil = flux_NOx_soil = flux_N2O_soil = flux_N2_soil = 0.0;
	flux_NH3_fire = flux_NOx_fire = flux_N2O_fire = flux_N2_fire = 0.0;
	flux_ntot = flux_nharvest = flux_nseed = 0.0;

	// Nitrogen in soil
	NH4_mass = NH3_mass = NO3_mass = 0.0;
	NO_mass = NO2_mass = N2O_mass = N2_mass = 0.0;
	NH3_mass_inc = NO_mass_inc = N2O_mass_inc = N2_mass_inc = 0.0;
	gross_nitrif = gross_denitrif = net_nitrif = net_denitrif= 0.0;

	double c_org_leach_lc[NLANDCOVERTYPES];

	for (int i=0; i<NLANDCOVERTYPES; i++) {
		c_org_leach_lc[i]=0.0;
	}

	// Sum C fluxes, dead C pools and runoff across patches

	Gridcell::iterator gc_itr = gridcell.begin();

	// Loop through Stands
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		stand.firstobj();

		//Loop through Patches
		while (stand.isobj) {
			Patch& patch = stand.getobj();

			double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

			flux_veg+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
			flux_repr+=-patch.fluxes.get_annual_flux(Fluxes::REPRC)*to_gridcell_average;
			flux_soil+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
			flux_fire+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
			flux_est+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
			flux_seed+=patch.fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
			flux_charvest+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;

			flux_nseed+=patch.fluxes.get_annual_flux(Fluxes::SEEDN)*to_gridcell_average;
			flux_nharvest+=patch.fluxes.get_annual_flux(Fluxes::HARVESTN)*to_gridcell_average;
			// N fluxes
			flux_NH3_fire+=patch.fluxes.get_annual_flux(Fluxes::NH3_FIRE)*to_gridcell_average;
			flux_NOx_fire+=patch.fluxes.get_annual_flux(Fluxes::NOx_FIRE)*to_gridcell_average;
			flux_N2O_fire+=patch.fluxes.get_annual_flux(Fluxes::N2O_FIRE)*to_gridcell_average;
			flux_N2_fire+=patch.fluxes.get_annual_flux(Fluxes::N2_FIRE)*to_gridcell_average;
			flux_NH3_soil+=patch.fluxes.get_annual_flux(Fluxes::NH3_SOIL)*to_gridcell_average;
			flux_NOx_soil+=patch.fluxes.get_annual_flux(Fluxes::NO_SOIL)*to_gridcell_average;
			flux_N2O_soil+=patch.fluxes.get_annual_flux(Fluxes::N2O_SOIL)*to_gridcell_average;
			flux_N2_soil+=patch.fluxes.get_annual_flux(Fluxes::N2_SOIL)*to_gridcell_average;	
			
			//Soil N flux from ntransform.cpp
			NH3_mass_inc+=patch.fluxes.get_annual_flux(Fluxes::NH3_SOIL)*to_gridcell_average;
			NO_mass_inc +=patch.fluxes.get_annual_flux(Fluxes::NO_SOIL)*to_gridcell_average;
			N2O_mass_inc+=patch.fluxes.get_annual_flux(Fluxes::N2O_SOIL)*to_gridcell_average;
			N2_mass_inc+=patch.fluxes.get_annual_flux(Fluxes::N2_SOIL)*to_gridcell_average;
			gross_nitrif+=patch.fluxes.get_annual_flux(Fluxes::GROSS_NITRIF)*to_gridcell_average;
			net_nitrif+=patch.fluxes.get_annual_flux(Fluxes::NET_NITRIF)*to_gridcell_average;
			gross_denitrif+=patch.fluxes.get_annual_flux(Fluxes::GROSS_DENITRIF)*to_gridcell_average;
			net_denitrif+=patch.fluxes.get_annual_flux(Fluxes::NET_DENITRIF)*to_gridcell_average;
			
			NH4_mass+=patch.soil.NH4_mass * to_gridcell_average;
			NO3_mass+=patch.soil.NO3_mass * to_gridcell_average;
			NO2_mass+=patch.soil.NO2_mass * to_gridcell_average;
			NO_mass+=patch.soil.NO_mass * to_gridcell_average;
			N2O_mass+=patch.soil.N2O_mass * to_gridcell_average;
			N2_mass+=patch.soil.N2_mass * to_gridcell_average;

			//Fire N flux
			flux_ntot+=(patch.fluxes.get_annual_flux(Fluxes::NH3_FIRE) +
						patch.fluxes.get_annual_flux(Fluxes::NOx_FIRE) +
						patch.fluxes.get_annual_flux(Fluxes::N2O_FIRE) +
						patch.fluxes.get_annual_flux(Fluxes::N2_FIRE) +
			//Soil N flux
					patch.fluxes.get_annual_flux(Fluxes::NH3_SOIL) +
					patch.fluxes.get_annual_flux(Fluxes::NO_SOIL) +
					patch.fluxes.get_annual_flux(Fluxes::N2O_SOIL) +
					patch.fluxes.get_annual_flux(Fluxes::N2_SOIL)) * to_gridcell_average;

			c_fast+=patch.soil.cpool_fast*to_gridcell_average;
			c_slow+=patch.soil.cpool_slow*to_gridcell_average;

			//Sum slow pools of harvested products
			if(run_landcover && ifslowharvestpool) {
				for (int q=0;q<npft;q++) {
					Patchpft& patchpft=patch.pft[q];
					c_harv_slow+=patchpft.harvested_products_slow*to_gridcell_average;
					n_harv_slow+=patchpft.harvested_products_slow_nmass*to_gridcell_average;
				}
			}

			surfrunoff_gridcell+=patch.asurfrunoff*to_gridcell_average;
			drainrunoff_gridcell+=patch.adrainrunoff*to_gridcell_average;
			baserunoff_gridcell+=patch.abaserunoff*to_gridcell_average;
			runoff_gridcell += patch.arunoff*to_gridcell_average;
			wetland_water_added_gridcell += patch.awetland_water_added*to_gridcell_average;

			// Fire return time
			if (!patch.has_fires() || patch.fireprob < 0.001) {
				firert_gridcell+=1000.0 * to_gridcell_average; // Set a limit of 1000 years
			}
			else {
				firert_gridcell+=(1.0/patch.fireprob) * to_gridcell_average;
				burned_area_gridcell+=patch.fireprob * to_gridcell_average;
			}

			aNH4dep_gridcell += gridcell.aNH4dep  * to_gridcell_average;
			aNO3dep_gridcell += gridcell.aNO3dep  * to_gridcell_average;
			anfert_gridcell += patch.anfert * to_gridcell_average;
			anmin_gridcell += patch.soil.anmin * to_gridcell_average;
			animm_gridcell += patch.soil.animmob * to_gridcell_average;
			anfix_gridcell += patch.soil.anfix * to_gridcell_average;
			n_min_leach_gridcell += patch.soil.aminleach * to_gridcell_average;
			n_org_leach_gridcell += patch.soil.aorgNleach * to_gridcell_average;
			c_org_leach_gridcell += patch.soil.aorgCleach * to_gridcell_average;
			availn += (patch.soil.NH4_mass + patch.soil.NO3_mass + patch.soil.snowpack_NH4_mass + patch.soil.snowpack_NO3_mass)
			           * to_gridcell_average;

			c_org_leach_lc[stand.landcover] += patch.soil.aorgCleach * to_gridcell_average;

			for (int r = 0; r < NSOMPOOL; r++) {

				if(r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
					surfsoillitterc += patch.soil.sompool[r].cmass * to_gridcell_average;
					surfsoillittern += patch.soil.sompool[r].nmass * to_gridcell_average;
				}
				else if (r == SURFFWD || r == SURFCWD) {
					cwdc += patch.soil.sompool[r].cmass * to_gridcell_average;
					cwdn += patch.soil.sompool[r].nmass * to_gridcell_average;
				}
				else {
					centuryc += patch.soil.sompool[r].cmass * to_gridcell_average;
					centuryn += patch.soil.sompool[r].nmass * to_gridcell_average;
				}
			}

			// Monthly output variables

			for (m=0;m<12;m++) {
				maet[m] += patch.maet[m]*to_gridcell_average;
				mpet[m] += patch.mpet[m]*to_gridcell_average;
				mevap[m] += patch.mevap[m]*to_gridcell_average;
				mintercep[m] += patch.mintercep[m]*to_gridcell_average;
				mrunoff[m] += patch.mrunoff[m]*to_gridcell_average;
				mrh[m] += patch.fluxes.get_monthly_flux(Fluxes::SOILC, m)*to_gridcell_average;
				mwcont_upper[m] += patch.soil.mwcont[m][0]*to_gridcell_average;
				mwcont_lower[m] += patch.soil.mwcont[m][1]*to_gridcell_average;

				mgpp[m] += patch.fluxes.get_monthly_flux(Fluxes::GPP, m)*to_gridcell_average;
				mra[m] += patch.fluxes.get_monthly_flux(Fluxes::RA, m)*to_gridcell_average;

				miso[m]+=patch.fluxes.get_monthly_flux(Fluxes::ISO, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_APIN, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_LIMO, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_TRIC, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_BPIN, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_MYRC, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_SABI, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_CAMP, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_TBOC, m)*to_gridcell_average;
				mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_OTHR, m)*to_gridcell_average;

				mmon_mt1[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_APIN, m)*to_gridcell_average;
				mmon_mt1[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_LIMO, m)*to_gridcell_average;
				mmon_mt1[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_TRIC, m)*to_gridcell_average;
				mmon_mt2[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_BPIN, m)*to_gridcell_average;
				mmon_mt2[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_MYRC, m)*to_gridcell_average;
				mmon_mt2[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_SABI, m)*to_gridcell_average;
				mmon_mt2[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_CAMP, m)*to_gridcell_average;
				mmon_mt2[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_TBOC, m)*to_gridcell_average;
				mmon_mt2[m]+=patch.fluxes.get_monthly_flux(Fluxes::MT_OTHR, m)*to_gridcell_average;
				
				msnowdepth[m] += patch.soil.msnowdepth[m]*M_PER_MM*to_gridcell_average;
				mald[m] += patch.soil.mthaw[m]*M_PER_MM*to_gridcell_average;			// mm to m
				mwtp[m] += patch.soil.mwtp[m]*M_PER_MM*to_gridcell_average;				// mm to m
				
				mch4[m] += patch.fluxes.get_monthly_flux(Fluxes::CH4C, m)*to_gridcell_average;				// g CH4-C/m2 - as CH4 is in gC, but CO2 fluxes are kgC
				mch4_diff[m] += patch.fluxes.get_monthly_flux(Fluxes::CH4C_DIFF, m)*to_gridcell_average;	// g CH4-C/m2
				mch4_plant[m] += patch.fluxes.get_monthly_flux(Fluxes::CH4C_PLAN, m)*to_gridcell_average;	// g CH4-C/m2
				mch4_ebull[m] += patch.fluxes.get_monthly_flux(Fluxes::CH4C_EBUL, m)*to_gridcell_average;	// g CH4-C/m2
				
				for (int sl = 0; sl < SOILTEMPOUT; sl++) {
					msoilt[m][sl] += patch.soil.T_soil_monthly[m][sl] * to_gridcell_average;
				}

			}

			maxald_gridcell += patch.soil.maxthawdepththisyear*M_PER_MM*to_gridcell_average;	// mm to m
			
			// Calculate monthly NPP and LAI

			Vegetation& vegetation = patch.vegetation;

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();

				// guess2008 - alive check added
				if (indiv.id != -1 && indiv.alive) {

					for (m=0;m<12;m++) {
						mlai[m] += indiv.mlai[m] * to_gridcell_average;
					}

				} // alive?

				vegetation.nextobj();

			} // while/vegetation loop
			stand.nextobj();
		} // patch loop
		++gc_itr;
	} // stand loop

	// In contrast to annual NEE, monthly NEE does not include fire
	// or establishment fluxes
	for (m=0;m<12;m++) {
		mnpp[m] = mgpp[m] - mra[m];
		mnee[m] = mrh[m] - mnpp[m];
	}

	// Print gridcell totals to files

	// Determine total leaf C:N ratio
	double cton_leaf_gridcell = limited_cton(cmass_leaf_gridcell, nmass_leaf_gridcell);

	// Determine total vmax nitrogen limitation
	if (cmass_leaf_gridcell > 0.0) {
		vmaxnlim_gridcell /= cmass_leaf_gridcell;
	}

	outlimit(out,out_cmass,					cmass_gridcell);
	outlimit(out,out_anpp,					anpp_gridcell);
	outlimit(out,out_agpp,					agpp_gridcell);
	outlimit(out,out_fpc,					fpc_gridcell);
	outlimit(out,out_aaet,					aaet_gridcell);
	outlimit(out,out_dens,					dens_gridcell);
	outlimit(out,out_lai,					lai_gridcell);
	outlimit(out,out_clitter,				clitter_gridcell);
	outlimit(out,out_aburned_area,				gridcell.annual_burned_area);
	outlimit(out,out_simfireanalysis,			gridcell.simfire_biome);
	outlimit(out,out_simfireanalysis,			gridcell.max_nesterov);
	outlimit(out,out_simfireanalysis,			gridcell.pop_density);
	outlimit(out,out_simfireanalysis,			gridcell.simfire_region);
	outlimit(out,out_firert,				firert_gridcell);
	outlimit(out,out_firert,				burned_area_gridcell);
	outlimit(out,out_runoff,				surfrunoff_gridcell);
	outlimit(out,out_runoff,				drainrunoff_gridcell);
	outlimit(out,out_runoff,				baserunoff_gridcell);
	outlimit(out,out_runoff,				runoff_gridcell);
	outlimit(out,out_wetland_water_added,	wetland_water_added_gridcell);
	
	outlimit(out,out_aiso,		aiso_gridcell);
	outlimit(out,out_amon,		amon_gridcell);
	outlimit(out,out_amon_mt1,  amon_mt1_gridcell);
	outlimit(out,out_amon_mt2,  amon_mt2_gridcell);

	outlimit(out,out_nmass,    (nmass_gridcell + nlitter_gridcell) * M2_PER_HA);
	outlimit(out,out_cton_leaf, cton_leaf_gridcell);
	outlimit(out,out_vmaxnlim,  vmaxnlim_gridcell);
	outlimit(out,out_nuptake,   nuptake_gridcell * M2_PER_HA);
	outlimit(out,out_nlitter,   nlitter_gridcell * M2_PER_HA);

	outlimit(out,out_nsources, aNH4dep_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, aNO3dep_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, anfix_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, anfert_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, (aNH4dep_gridcell + aNO3dep_gridcell + anfix_gridcell + anfert_gridcell) * M2_PER_HA);
	outlimit(out,out_nsources, anmin_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, animm_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, (anmin_gridcell - animm_gridcell) * M2_PER_HA);
	outlimit(out,out_nsources, (anmin_gridcell - animm_gridcell + aNH4dep_gridcell + aNO3dep_gridcell +
				anfix_gridcell + anfert_gridcell) * M2_PER_HA);

	// Print landcover totals to files
	if (run_landcover) {
		for(int i=0;i<NLANDCOVERTYPES;i++) {
			if(run[i]) {
				outlimit(out,out_cmass, landcover_cmass[i]);
				outlimit(out,out_anpp,  landcover_anpp[i]);
				outlimit(out,out_agpp,  landcover_agpp[i]);
				outlimit(out,out_fpc,   landcover_fpc[i]);
				outlimit(out,out_aaet,  landcover_aaet[i]);
				outlimit(out,out_dens,  landcover_densindiv_total[i]);
				outlimit(out,out_lai,   landcover_lai[i]);
				outlimit(out,out_clitter, landcover_clitter[i]);
				outlimit(out,out_aiso,  landcover_aiso[i]);
				outlimit(out,out_amon,  landcover_amon[i]);
				outlimit(out,out_amon_mt1,  landcover_amon_mt1[i]);
				outlimit(out,out_amon_mt2,  landcover_amon_mt2[i]);

				double landcover_cton_leaf = limited_cton(landcover_cmass_leaf[i],
						landcover_nmass_leaf[i]);

				if (landcover_cmass_leaf[i] > 0.0) {
					landcover_vmaxnlim[i] /= landcover_cmass_leaf[i];
				}

				outlimit(out,out_nmass,     (landcover_nmass[i] + landcover_nlitter[i]) * M2_PER_HA);
				outlimit(out,out_cton_leaf, landcover_cton_leaf);
				outlimit(out,out_vmaxnlim,  landcover_vmaxnlim[i]);
				outlimit(out,out_nuptake,   landcover_nuptake[i] * M2_PER_HA);
				outlimit(out,out_nlitter,   landcover_nlitter[i] * M2_PER_HA);
			}
		}
	}

	// Print monthly output variables
	for (m=0;m<12;m++) {
		outlimit(out,out_mnpp,         mnpp[m]);
		outlimit(out,out_mlai,         mlai[m]);
		outlimit(out,out_mgpp,         mgpp[m]);
		outlimit(out,out_mra,          mra[m]);
		outlimit(out,out_maet,         maet[m]);
		outlimit(out,out_mpet,         mpet[m]);
		outlimit(out,out_mevap,        mevap[m]);
		outlimit(out,out_mrunoff,      mrunoff[m]);
		outlimit(out,out_mintercep,    mintercep[m]);
		outlimit(out,out_mrh,          mrh[m]);
		outlimit(out,out_mnee,         mnee[m]);
		outlimit(out,out_mwcont_upper, mwcont_upper[m]);
		outlimit(out,out_mwcont_lower, mwcont_lower[m]);
		outlimit(out,out_miso,         miso[m]);
		outlimit(out,out_mmon,         mmon[m]);
		outlimit(out,out_mmon_mt1,     mmon_mt1[m]);
		outlimit(out,out_mmon_mt2,     mmon_mt2[m]);
		outlimit(out,out_mburned_area, (float)gridcell.monthly_burned_area[m]);

		aaet += maet[m];
		apet += mpet[m];
		aevap += mevap[m];
		arunoff += mrunoff[m];
		aintercep += mintercep[m];

		// Arctic and wetland output
		const int layer_ix_25cm = 2; // Layer index for 25cm soil depth. It could depend on the thickness of the layers in future updates.  
		outlimit(out,out_msoiltempdepth5, msoilt[m][0]);
		outlimit(out,out_msoiltempdepth15, msoilt[m][1]);
		outlimit(out,out_msoiltempdepth25, msoilt[m][layer_ix_25cm]);
		outlimit(out,out_msoiltempdepth35, msoilt[m][3]);
		outlimit(out,out_msoiltempdepth45, msoilt[m][4]);
		outlimit(out,out_msoiltempdepth55, msoilt[m][5]);
		outlimit(out,out_msoiltempdepth65, msoilt[m][6]);
		outlimit(out,out_msoiltempdepth75, msoilt[m][7]);
		outlimit(out,out_msoiltempdepth85, msoilt[m][8]);
		outlimit(out,out_msoiltempdepth95, msoilt[m][9]);
		outlimit(out,out_msoiltempdepth105, msoilt[m][10]);
		outlimit(out,out_msoiltempdepth115, msoilt[m][11]);
		outlimit(out,out_msoiltempdepth125, msoilt[m][12]);
		outlimit(out,out_msoiltempdepth135, msoilt[m][13]);
		outlimit(out,out_msoiltempdepth145, msoilt[m][14]);
		outlimit(out,out_mch4, mch4[m]);
		outlimit(out,out_mch4diff, mch4_diff[m]);
		outlimit(out,out_mch4plan, mch4_plant[m]);
		outlimit(out,out_mch4ebull, mch4_ebull[m]);		
		outlimit(out,out_msnow, msnowdepth[m]);
		outlimit(out,out_mwtp, mwtp[m]);
		outlimit(out,out_mald, mald[m]);
	}

	outlimit(out,out_mald, maxald_gridcell); // [m]

	// Graphical output every PLOT_INTERVAL years
	// (Windows shell only - no effect otherwise)

	if (!(date.year%PLOT_INTERVAL)) {
		if(gridcell.nbr_stands() > 0)	//Fixed bug here if no stands were present.
		{
			Stand& stand = gridcell[0];
			plot("C flux [kgC/m2/yr]","veg",  date.year, flux_veg);
			plot("C flux [kgC/m2/yr]","repr", date.year, flux_repr);
			plot("C flux [kgC/m2/yr]","soil", date.year, flux_soil);
			plot("C flux [kgC/m2/yr]","fire", date.year, flux_fire);
			plot("C flux [kgC/m2/yr]","est",  date.year, flux_est);
			plot("C flux [kgC/m2/yr]","NEE",  date.year, flux_veg + flux_repr + flux_soil + flux_fire + flux_est);

			if (!ifcentury) {
				plot("Soil C [kgC/m2]","slow", date.year, stand[0].soil.cpool_slow);
				plot("Soil C [kgC/m2]","fast", date.year, stand[0].soil.cpool_fast);
			}
			else {
				plot("N flux [kgN/ha/yr]","fix",   date.year, -anfix_gridcell * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","dep",   date.year, -(aNH4dep_gridcell - aNO3dep_gridcell) * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","fert",  date.year, -anfert_gridcell * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","leach", date.year, (n_min_leach_gridcell + n_org_leach_gridcell) * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","emissions",  date.year, flux_ntot * M2_PER_HA);

				plot("N flux [kgN/ha/yr]","NEE",   date.year, (flux_ntot + n_min_leach_gridcell + n_org_leach_gridcell -
					(anfix_gridcell + aNH4dep_gridcell + aNO3dep_gridcell + anfert_gridcell)) * M2_PER_HA);

				plot("N mineralization [kgN/ha/yr]","N", date.year, (anmin_gridcell - animm_gridcell) * M2_PER_HA);

				plot("Soil C [kgC/m2]","fine litter",   date.year, surfsoillitterc);
				plot("Soil C [kgC/m2]","coarse litter", date.year, cwdc);
				plot("Soil C [kgC/m2]","soil",          date.year, centuryc);
				plot("Soil C [kgC/m2]","total",         date.year, surfsoillitterc + cwdc + centuryc);

				plot("Soil N [kgN/m2]","fine litter",   date.year, surfsoillittern);
				plot("Soil N [kgN/m2]","coarse litter", date.year, cwdn);
				plot("Soil N [kgN/m2]","soil",          date.year, centuryn);
				plot("Soil N [kgN/m2]","total",         date.year, surfsoillittern + cwdn + centuryn);
			}

			plot("H2O flux [mm/yr]", "transp", date.year, aaet);
			plot("H2O flux [mm/yr]", "intercep", date.year, aintercep);
			plot("H2O flux [mm/yr]", "evap", date.year, aevap);
			plot("H2O flux [mm/yr]", "runoff", date.year, arunoff);
			plot("H2O flux [mm/yr]", "PET", date.year, apet);
		}
	}

    // Write fluxes to file

	Landcover& lc = gridcell.landcover;

	outlimit(out,out_cflux, flux_veg);
	outlimit(out,out_cflux, -flux_repr);
	outlimit(out,out_cflux, flux_soil + c_org_leach_gridcell);
	outlimit(out,out_cflux, flux_fire);
	outlimit(out,out_cflux, flux_est);
	if (run_landcover) {
			outlimit(out,out_cflux, flux_seed);
			outlimit(out,out_cflux, flux_charvest);
			outlimit(out,out_cflux, lc.acflux_landuse_change);
			outlimit(out,out_cflux, lc.acflux_harvest_slow);
	}
	outlimit(out,out_cflux, flux_veg - flux_repr + flux_soil + flux_fire + flux_est + c_org_leach_gridcell +
			flux_seed + flux_charvest + lc.acflux_landuse_change + lc.acflux_harvest_slow);

	outlimit(out,out_doc, (c_org_leach_gridcell) * M2_PER_HA);

	if (run_landcover) {
		for(int i=0;i<NLANDCOVERTYPES;i++) {
			if(run[i]) {
				outlimit(out,out_doc, c_org_leach_lc[i] * M2_PER_HA);
			}
		}
	}


	outlimit(out,out_nflux, -aNH4dep_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, -aNO3dep_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, -anfix_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, -anfert_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, flux_ntot * M2_PER_HA);
	outlimit(out,out_nflux, (n_min_leach_gridcell + n_org_leach_gridcell) * M2_PER_HA);

	if (run_landcover) {
			outlimit(out,out_nflux, flux_nseed * M2_PER_HA);
			outlimit(out,out_nflux, flux_nharvest * M2_PER_HA);
			outlimit(out,out_nflux, lc.anflux_landuse_change * M2_PER_HA);
			outlimit(out,out_nflux, lc.anflux_harvest_slow * M2_PER_HA);
	}
	outlimit(out,out_nflux, (flux_nharvest + lc.anflux_landuse_change +
				lc.anflux_harvest_slow + flux_nseed + flux_ntot +
				n_min_leach_gridcell + n_org_leach_gridcell -
					 (aNH4dep_gridcell + aNO3dep_gridcell + anfix_gridcell + anfert_gridcell)) * M2_PER_HA);

	// CPOOL Write cpool to file

	outlimit(out,out_cpool, cmass_gridcell);
	if (!ifcentury) {
		outlimit(out,out_cpool, clitter_gridcell);
		outlimit(out,out_cpool, c_fast);
		outlimit(out,out_cpool, c_slow);
	}
	else {
		outlimit(out,out_cpool, clitter_gridcell + surfsoillitterc + cwdc);
		outlimit(out,out_cpool, centuryc);
	}

	if (run_landcover && ifslowharvestpool) {
		outlimit(out,out_cpool, c_harv_slow);
	}

	// Calculate total cpool, starting with cmass and litter...
	double cpool_total = cmass_gridcell + clitter_gridcell;

	// Add SOM pools
	if (!ifcentury) {
		cpool_total += c_fast + c_slow;
	}
	else {
		cpool_total += centuryc + surfsoillitterc + cwdc;
	}

	// Add slow harvest pool if needed
	if (run_landcover && ifslowharvestpool) {
		cpool_total += c_harv_slow;
	}

	outlimit(out,out_cpool, cpool_total);

	// NPOOL Write npool to file
	
	if (ifcentury) {
		outlimit(out,out_npool, nmass_gridcell + nlitter_gridcell);
		outlimit(out,out_npool, surfsoillittern + cwdn);
		outlimit(out,out_npool, centuryn + availn);

		if(run_landcover && ifslowharvestpool) {
			outlimit(out,out_npool, n_harv_slow);
			outlimit(out,out_npool, (nmass_gridcell + nlitter_gridcell + surfsoillittern + cwdn + centuryn + availn + n_harv_slow));
		}
		else {
			outlimit(out,out_npool, (nmass_gridcell + nlitter_gridcell + surfsoillittern + cwdn + centuryn + availn));
		}
	}

	outlimit(out,out_ngases, flux_NH3_fire   * M2_PER_HA);
	outlimit(out,out_ngases, flux_NH3_soil   * M2_PER_HA);
	outlimit(out,out_ngases, flux_NOx_fire   * M2_PER_HA);
	outlimit(out,out_ngases, flux_NOx_soil   * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2O_fire   * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2O_soil   * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2_fire    * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2_soil    * M2_PER_HA);
	outlimit(out,out_ngases, flux_ntot  * M2_PER_HA);
	// Volatile nitrogen from soil
	outlimit(out,out_soil_nflux, NH3_mass_inc * M2_PER_HA);
	outlimit(out,out_soil_nflux, NO_mass_inc * M2_PER_HA);
	outlimit(out,out_soil_nflux, N2O_mass_inc * M2_PER_HA);
	outlimit(out,out_soil_nflux, N2_mass_inc * M2_PER_HA);
	// Nitrification/denitrification fluxes
	//outlimit(out,out_soil_nflux, gross_nitrif * M2_PER_HA);
	//outlimit(out,out_soil_nflux, net_nitrif * M2_PER_HA);
	//outlimit(out,out_soil_nflux, gross_denitrif * M2_PER_HA);
	//outlimit(out,out_soil_nflux, net_denitrif * M2_PER_HA);
		
	// Ackumulated nitrogen in pools		
	outlimit(out,out_soil_npool, NH4_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, NO3_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, NO2_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, NO_mass  * M2_PER_HA);
	outlimit(out,out_soil_npool, N2O_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, N2_mass  * M2_PER_HA);


	// Output of tree stand age structure, monthly soil water and 3D vegetation view
	// (Windows shell only - no effect otherwise)

	if (vegmode==COHORT || vegmode==INDIVIDUAL) {

		if (!(date.year%PLOT_UPDATE_INTERVAL)) {

			double* densindiv=NULL;
			int nageclass;
			get_stand_age_structure(gridcell, densindiv, nageclass, false);

			if (nageclass) {

				densindiv = new double[(npft+1)*nageclass];
				if (densindiv) {

					resetwindow("Age structure [indiv/ha]");

					get_stand_age_structure(gridcell, densindiv, nageclass, true);

					pftlist.firstobj();
					while (pftlist.isobj) {
						Pft& pft = pftlist.getobj();

						if (pft.lifeform == TREE) {

							for (c = 0; c<nageclass; c++)
								plot("Age structure [indiv/ha]", pft.name,
								c * estinterval + estinterval*0.5,
								densindiv[pft.id*nageclass+c]*1e4); // includes conversion from /m2 --> /ha
						}

						pftlist.nextobj();
					}

					delete[] densindiv;
				}
			}
		}
	}

	if (!(date.year % PLOT_UPDATE_INTERVAL)) {

		resetwindow("Soil water [%AWC]");

		for (m = 0; m < 12; m++) {
			plot("Soil water [%AWC]", "upper", m, mwcont_upper[m]);
			plot("Soil water [%AWC]", "lower", m, mwcont_lower[m]);
		}
	}

	if (vegmode == COHORT) {
		if (!(date.year%VEG3D_UPDATE_INTERVAL)) {
			if (date.year == 0) open3d();
			output_vegetation(gridcell, pftlist);
		}
	}

}

/// Output of simulation results at the end of each day
/** This function does not have to provide any information to the framework.
  */
void CommonOutput::outdaily(Gridcell& gridcell) {
}

} // namespace
