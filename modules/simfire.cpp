///////////////////////////////////////////////////////////////////////////////////////
/// \file simfire.cpp
/// \brief SIMFIRE burned area simulation by W. Knorr
///
/// \author Lars Nieradzik
/// $Date: 2017-01-24 16:02:51 +0100 (Tue, 24 Jan 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"
#include "simfire.h"
#include "simfire_input.h"

const int NFIREBIOMES = 9;

// Annually update SIMFIRE BIOME for a patch from current vegetation
/* Computes current SIMFIRE biome for this
 * gridcell depending on the last <n_year_biomeavg> years of
 * vegetation.
 */
int update_fire_biome(Patch& patch, double lat) {
	
	/* SIMFIRE BIOMES:
	 * 0 no veg/no data
	 * 1 Cropland/Urban/Natural Vegetation Mosaic (IGBP 12-14)
	 * 2 Needleleaf forest (IGBP 1,3): >60% cover, height>2m
	 * 3 Broadleaf forest (IGBP 2,4): >60% cover, height>2m
	 * 4 Mixed forest (IGBP 4): >60% cover, height>2m, none >60%
	 * 5 Shrubland (IGBP 6,7 and latitude<50): >10% woody cover, height<2m
	 * 6 Savanna or Grassland (IGBP 8-10): herbaceous component present, <60% tree cover
	 * 7 Tundra (IGBP 6,7,16 and latitude>=50): height<2m
	 * 8 Barren or Sparsely Vegetated (IGBP 16 and latitude<50): <10% vegetation cover
	 */
	
	double fpar_grass      = 0.0; // grass fraction of all vegetation
	double fpar_needleleaf = 0.0; // fraction of needle-leaf tress
	double fpar_broadleaf  = 0.0; // fraction of broad-leaf trees
	double fpar_trop_broadleaf_raingreen  = 0.0; // fraction of TrBR (for african savannah)
	double fpar_shrubs     = 0.0; // fraction of woody vegetation that is shrubs
	double fpar_total      = 0.0; // total FPAR of all individuals
	int biome = 0;   // biome number

	// Obtain reference to Vegetation object for this patch
	Vegetation& vegetation=patch.vegetation;

	// Loop through individuals of this patch
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();
	
		if (indiv.id!=-1 && indiv.alive) {
	
			if (indiv.pft.lifeform==GRASS) {
				fpar_grass += indiv.fpar_leafon;
			}
			else { // tree or shrub
				if (indiv.height>=2.0) { // tree
					if (indiv.pft.leafphysiognomy==NEEDLELEAF) {
						fpar_needleleaf += indiv.fpar_leafon;
					}
					else { // broadleaf tree
						// For savannas this is accounted for separately
						if (indiv.pft.phenology==RAINGREEN){
							fpar_trop_broadleaf_raingreen+=indiv.fpar_leafon; 
						}
						fpar_broadleaf += indiv.fpar_leafon;
					}
				}
				else { // shrub
					fpar_shrubs += indiv.fpar_leafon;
				}
			}
			fpar_total += indiv.fpar_leafon;
		}
		vegetation.nextobj(); // ... on to next individual
	}
	
	if ( fpar_total < 0.00000001 ) {
		return SF_NOVEG; // no data
	}

	// Re-normalize
	fpar_grass                     /= fpar_total;
	fpar_needleleaf                /= fpar_total;
	fpar_broadleaf                 /= fpar_total;
	fpar_trop_broadleaf_raingreen  /= fpar_total;
	fpar_shrubs                    /= (1.00001-fpar_grass);

	// Save current 
	int idx = date.year % N_YEAR_BIOMEAVG;
	patch.fapar_total_avg[idx] = fpar_total;
	patch.fapar_grass_avg[idx] = fpar_grass;
	patch.fapar_ndlt_avg [idx] = fpar_needleleaf;
	patch.fapar_brlt_avg [idx] = fpar_broadleaf;
	patch.fapar_trbr_avg [idx] = fpar_trop_broadleaf_raingreen;
	patch.fapar_shrub_avg[idx] = fpar_shrubs;
	
	// Generate running avereage 
	fpar_total      = 0.;
	fpar_grass      = 0.;
	fpar_needleleaf = 0.;
	fpar_broadleaf  = 0.;
	fpar_trop_broadleaf_raingreen = 0.;
	fpar_shrubs     = 0.;
	for (int i = 0; i<N_YEAR_BIOMEAVG; i++) {
		fpar_total      += patch.fapar_total_avg[i];
		fpar_grass      += patch.fapar_grass_avg[i];
		fpar_needleleaf += patch.fapar_ndlt_avg [i];
		fpar_broadleaf  += patch.fapar_brlt_avg [i];
		fpar_trop_broadleaf_raingreen  += patch.fapar_trbr_avg [i];
		fpar_shrubs     += patch.fapar_shrub_avg[i];
	}
	fpar_total      /=  (double)N_YEAR_BIOMEAVG;
	fpar_grass      /=  (double)N_YEAR_BIOMEAVG;
	fpar_needleleaf /=  (double)N_YEAR_BIOMEAVG;
	fpar_broadleaf  /=  (double)N_YEAR_BIOMEAVG;
	fpar_trop_broadleaf_raingreen  /=  (double)N_YEAR_BIOMEAVG;
	fpar_shrubs     /=  (double)N_YEAR_BIOMEAVG;

	// Set Fire-Biome
	if (fpar_total < 0.5 && fabs(lat) < 50.0) {
		biome = SF_BARREN;
	} 
	else if (fpar_total < 0.35   && fabs(lat) >= 50.0) {
		biome = SF_TUNDRA;
	} 
	else if (patch.stand.landcover == CROPLAND) {
		biome = SF_CROP;
	} 
	else if (fpar_shrubs >= 0.9 && fabs(lat) < 50.0) {
		biome = SF_SHRUBS;
	} 
	else if (fpar_shrubs >= 0.4 && fabs(lat) >= 50.0) {
		biome = SF_TUNDRA;
	} 
	else if (fpar_needleleaf >= 0.8) {
		biome = SF_NEEDLELEAF;
	} 
	else if (fpar_broadleaf >= 0.8) {
		if (fpar_trop_broadleaf_raingreen >=0.3) {
			biome = SF_SAVANNA;
		}
		else {
			biome = SF_BROADLEAF;
		}
	}
	else if (fpar_grass >= 0.3) {
		biome = SF_SAVANNA;
	} 
	else if (fpar_broadleaf > 0.2 && fpar_needleleaf > 0.2){
		biome = SF_MIXED_FOREST;
	}
	else {
		biome = SF_SAVANNA;
	}

	return biome;
}

// Compute dominant SIMFIRE biome for gridcell
/* Computes current SIMFIRE biome for this
 * gridcell depending on the last <n_year_biomeavg> years of
 * vegetation.
 */
void simfire_biome_mapping(Gridcell& gridcell) {

	std::vector<int> biomes;
	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		
		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			if (patch.stand.landcover != BARREN && patch.stand.landcover != CROPLAND && patch.stand.landcover != URBAN) {
				biomes.push_back(update_fire_biome(patch, gridcell.get_lat()));
			}
			stand.nextobj();
		}
		++gc_itr;
	}
		
	int count[NFIREBIOMES];
	int biome;
	int count_max=0; // maximum of 'count'
	
	// Find and save most common biome number
	for (biome=0; biome < NFIREBIOMES; biome++) {
		count[biome] = 0;
	}
	for (int idx = 0; idx < (int)biomes.size(); idx++) {
		count[biomes[idx]]++;
	}
	for (biome=0; biome < NFIREBIOMES; biome++) {
		count_max = max(count_max, count[biome]);
	}
	
	int biome_index = 0;
	for (biome = 0; biome < NFIREBIOMES && count[biome] < count_max; biome++) {
		biome_index++;
	}

	// Determine SIMFIRE-Biome for whole gridcell
	double frac_nat = gridcell.landcover.frac[NATURAL] + gridcell.landcover.frac[PEATLAND] + gridcell.landcover.frac[FOREST];
	double frac_barren = gridcell.landcover.frac[URBAN] + gridcell.landcover.frac[BARREN];
	if (frac_nat > 0.5) {
		gridcell.simfire_biome = biome_index;
	}
	else if (frac_barren > 0.5) {
		gridcell.simfire_biome = SF_BARREN;
	}
	else if (gridcell.landcover.frac[CROPLAND] > 0.5) {
		gridcell.simfire_biome = SF_NOVEG;
	}
	else if (gridcell.landcover.frac[CROPLAND] > 0.3) {
		gridcell.simfire_biome = SF_CROP;
	}
	else if ( gridcell.landcover.frac[CROPLAND] < 0.1 && frac_nat < 0.4) {
		gridcell.simfire_biome = SF_SHRUBS;
	}
	else {
		gridcell.simfire_biome = SF_SAVANNA;
	}
}

// Read SIMFIRE related data for a gridcell at beginning of simulation
/* Reads SIMFIRE relevant info from simfire_input.bin:
 * - Hyde 3.1 population density
 * - Monthly fire climatology
 */
void getsimfiredata(Gridcell& gridcell) {

	// Paths to SIMFIRE binaries
	xtring file_simfire = param["file_simfire"].str;
	
	// Open file, fill podp, monthly_burned_area and igbp_class for a gridcell
	// Fill static arrays/variables here
	SimfireInputArchive ark;
	
	if (!ark.open(file_simfire)) {
		fail("Could not open %s for input \n", (char*)file_simfire);
	}
	
	SimfireInput rec;

	// Make sure gridcell lat/lon is centered to LPJ-GUESS gridcell
	rec.lon = floor(gridcell.get_lon() * 2.) / 2. + 0.25;
	rec.lat = floor(gridcell.get_lat() * 2.) / 2. + 0.25;
	
	if (!ark.getindex(rec)) {
		ark.close();
		fail("Grid cell not found in %s \n", (char*)file_simfire);
	}

	// Found the record, get the values
	
	// Convert IGBP into simfire internal biomes
	simfire_biome_mapping(gridcell);
	
	// Monthly fire risk (W.Knorr)
	for (int m = 0; m < 12; m++) {
		gridcell.monthly_fire_risk[m] = rec.monthly_burned_area[m];
	}

	// Population density from HYDE 3.1
	for (int t=0; t<57; t++) {
		gridcell.hyde31_pop_density[t] = rec.pop_density[t];
	}		
	
	ark.close();
}

// Get this year's human population density
/* Computes population density from the Hyde 3.1 dataset.
 * Annual data is computed by linearly interpolating between the existing values.
 * Before 10000 BC the 10000 BC value is used, after 2005 linear extrapolation
 * using the change between the last two values is performed
 */
void simfire_update_pop_density(Gridcell& gridcell) {

	// Number of entries for Population data
	const int NPOPENTRIES = 57;

	// Years at which population-data is available in HYDE3.1
	const int POPTIME[NPOPENTRIES]  = {-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,
		-2000,-1000,0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,
		1500,1600,1700,1710,1720,1730,1740,1750,1760,1780,1790,1810,1820,1830,1840,
		1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,
		2000,2005};
	
	// Get calendar-year
	int cyear = date.get_calendar_year();

	// Find start and end year index of pop interpolation
	int idx = 0 ;    
	while (POPTIME[idx] < cyear && idx < NPOPENTRIES-1) {
		idx++;
	}

	double popd;
	if ( cyear <= POPTIME[0] ) {
		// Use first year's value (10000 BC) for earlier years.
		popd = gridcell.hyde31_pop_density[0];
	}
	else if ( cyear >= POPTIME[NPOPENTRIES-1] ) {
		// linearly extrapolate latest growth/decline
		popd = gridcell.hyde31_pop_density[NPOPENTRIES-1] + 
			(gridcell.hyde31_pop_density[NPOPENTRIES-1]-gridcell.hyde31_pop_density[NPOPENTRIES-2]) /
			(double)(POPTIME[NPOPENTRIES-1] - POPTIME[NPOPENTRIES-2]) * (double)(cyear-POPTIME[NPOPENTRIES-1]);
	}
	else {

		// Interpolate between two entries
		double interpf = (double)(cyear-POPTIME[idx-1]) /
			(double)( POPTIME[idx]-POPTIME[idx-1] );
		popd = (1. - interpf) * gridcell.hyde31_pop_density[idx-1] + 
			interpf * gridcell.hyde31_pop_density[idx];
	}

	gridcell.pop_density = max(0.,popd);
}

// Daily book-keeping for SIMFIRE-relevant variables
/* Updates SIMFIRE's Maximum Annual Mesterov Index
 * and running mean of max annual FPAR (from canexch.cpp).
 * Updates fire biome at beginning of the year.
 */
void simfire_accounting_gridcell(Gridcell& gridcell) {
	
	Climate& climate = gridcell.climate;
	
	// Absolute upper boundary for the accumulative nesterov index
	const double MAXIMUM_NESTEROV = 250000; 

	// Check if it is first day of simulation
	bool is_first_day = ( date.day == 0 && ( date.year == 0 || 
			       ( restart && date.year == state_year ) ) );

	if ( is_first_day ) {
		// Read monthly climatology and Hyde population-data from file 
		getsimfiredata(gridcell);
	}

	if (date.day == 0 ) {
		
		double avg = 0.;
		for(int i=0;i<AVG_INTERVAL_FAPAR;i++) {
			avg += gridcell.fapar_recent_max[i];
		}
		gridcell.fapar_annual_max = avg / (double) AVG_INTERVAL_FAPAR;

		// Finally (re)set this years max fapar
		gridcell.fapar_cur_max = 0.0;

		// Set global simfire region as fixed: Global=0
		gridcell.simfire_region = 0;

		// Determine SIMFIRE biome for this year
		simfire_biome_mapping(gridcell);

		// Update population density
		simfire_update_pop_density(gridcell);
	}
	// Multi-year accounting of maximum annual fapar	
	else if ( date.islastday && date.islastmonth ) {
		int a = date.year % AVG_INTERVAL_FAPAR;
		gridcell.fapar_recent_max[a] = gridcell.fapar_cur_max;
	}

	// Update running Maximum Nesterov index array at beginning of month  
	if ( date.dayofmonth == 0 ) {
		double mnest = 0.;
		for ( int i=0; i < 12; i++) 
			if ( gridcell.nesterov_monthly_max[i] > mnest ) {
				mnest = gridcell.nesterov_monthly_max[i];
			}
		gridcell.nesterov_max = mnest;
		gridcell.nesterov_monthly_max[date.month] = 0. ;
	}

	// Update current month's Maximum Nesterov index
	if (  gridcell.nesterov_monthly_max[date.month] < gridcell.nesterov_cur ) {
		gridcell.nesterov_monthly_max[date.month] = gridcell.nesterov_cur;
	}
	
	// Loop over patches for fpar
	int cnt= 0;
	double fapar_run_mean = 0.;
	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			if ( ! ( date.year == 0 && date.day == 0 ) ) {
				fapar_run_mean += (1. - patch.fpar_ff);
			}
			cnt += 1;
			stand.nextobj();
		}
		++gc_itr;
	}
	// Average over each patch
	fapar_run_mean /= (double) cnt;

	// Update the this years maximum
	gridcell.fapar_cur_max = max(fapar_run_mean, gridcell.fapar_cur_max);

	// Compute running Nesterov index
	if ( climate.prec >= 3. || climate.tmax - climate.tmin < 4. ) {
		gridcell.nesterov_cur = 0.0;
	}
	else {
		gridcell.nesterov_cur += ( climate.tmax - climate.tmin + 4. ) * climate.tmax ;
	}
	gridcell.nesterov_cur = min(gridcell.nesterov_cur,MAXIMUM_NESTEROV) ;

	// Finally update Max Annual Mesterov Index
	if (gridcell.nesterov_cur > gridcell.nesterov_max ) {
		gridcell.nesterov_max = gridcell.nesterov_cur ;
	}
}

// Calculate burned area in ha following Knorr 2014.
double simfire_burned_area(Gridcell& gridcell) {

	// Globally trained parameters 
	//                    Crop    NL     BL      Mixed  Shrub  Savana Tundra Barren  
	const double A[8] = { 0.110,  0.095, 0.092  ,0.127, 0.470, 0.889, 0.059, 0.113  }; 
	const double B = 0.905;  
	const double C = 0.860; 
	const double E = -0.0168; 
	const double SCALAR = 1.0e-5;

	// Return if improper biome-type
	if (gridcell.simfire_biome == SF_NOVEG) {
		return 0.;
	}

	// fPAR correction Knorr for use with LPJ-GUESS only
	const double FPAR_CORR1 = 0.428;
	const double FPAR_CORR2 = 0.148;
	double fpar_cor = FPAR_CORR1 * gridcell.fapar_annual_max + FPAR_CORR2 * gridcell.fapar_annual_max *
	  gridcell.fapar_annual_max;

	// Compute annual burned area
	double burned_area = A[gridcell.simfire_biome-1] *
		pow(fpar_cor, B) *
		pow((SCALAR * gridcell.nesterov_max), C) *
		exp(E * gridcell.pop_density);

	// Compute daily burned_area
	burned_area *= gridcell.monthly_fire_risk[date.month] /
		(double)date.ndaymonth[date.month];

	// For annual simfire firebrobability
	gridcell.simfire_annual_burned_area += burned_area;
	
	return burned_area;
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Knorr, W. et al., Impact of human population density on fire frequency at the 
//   global scale, BIOGEOSCIENCES, 11, 4, 2014, DOI: 10.5194/bg-11-1085-2014
