///////////////////////////////////////////////////////////////////////////////////////
/// \file soil.cpp
/// \brief LPJ-GUESS Soil class and related definitions
/// Implementation of member functions of class Soil.
/// The class Soil and its member functions and variables are declared in guess.h
///
/// \author Paul Miller
/// $Date: 2015-12-22 15:22:04 +0100 (Tue, 22 Dec 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"
#include "soil.h"
#include "parameters.h"


////////////////////////////////////////////////////////////////////////////////
// Implementation of Soil member functions
////////////////////////////////////////////////////////////////////////////////

void Soil::init_states() {

	// Initialises certain member variables

	alag = 0.0;
	exp_alag = 1.0;
	cpool_slow = 0.0;
	cpool_fast = 0.0;
	decomp_litter_mean = 0.0;
	k_soilfast_mean = 0.0;
	k_soilslow_mean = 0.0;
	wcont_evap = 0.0;
	snowpack = 0.0;
	orgleachfrac = 0.0;

	// Extra initialisation
	aorgCleach = 0.0;
	aorgNleach = 0.0;
	anfix = 0.0;
	aminleach = 0.0;

	lKorg = lKpeat = lKmin = lKwater = lKice = lKair = 0.0;

	mwcontupper = 0.0;
	mwcontlower = 0.0;

	for (int mth=0; mth<12; mth++) {
		mwcont[mth][0] = 0.0;
		mwcont[mth][1] = 0.0;
		fnuptake_mean[mth] = 0.0;
		morgleach_mean[mth] = 0.0;
		mminleach_mean[mth] = 0.0;
	}

	std::fill_n(dwcontupper, Date::MAX_YEAR_LENGTH, 0.0);
	std::fill_n(dwcontlower, Date::MAX_YEAR_LENGTH, 0.0);
	// for fire
	std::fill_n(dthaw, Date::MAX_YEAR_LENGTH, 0.0);


	/////////////////////////////////////////////////////
	// Initialise CENTURY pools

	// Set initial CENTURY pool N:C ratios
	// Parton et al 1993, Fig 4

	sompool[SOILMICRO].ntoc = 1.0 / 15.0;
	sompool[SURFHUMUS].ntoc = 1.0 / 15.0;
	sompool[SLOWSOM].ntoc = 1.0 / 20.0;
	sompool[SURFMICRO].ntoc = 1.0 / 20.0;

	// passive has a fixed value
	sompool[PASSIVESOM].ntoc = 1.0 / 9.0;
	
	NO2_mass = 0.0;
	NO2_mass_w = 0.0;
	NO2_mass_d = 0.0;
	NO_mass = 0.0;
	NO_mass_w = 0.0;
	NO_mass_d = 0.0;
	N2O_mass = 0.0;
	N2O_mass_w = 0.0;
	N2O_mass_d = 0.0;
	N2_mass = 0.0;
	NH4_mass = 0.0;
	NO3_mass = 0.0;
	NH4_input = 0.0;
	NO3_input = 0.0;
	anmin = 0.0;
	animmob = 0.0;
	aminleach = 0.0;
	aorgNleach = 0.0;
	aorgCleach = 0.0;
	anfix = 0.0;
	anfix_calc = 0.0;
	anfix_mean = 0.0;
	snowpack_NH4_mass = 0.0;
	snowpack_NO3_mass = 0.0;
	labile_carbon = 0.0;
	labile_carbon_w = 0.0;
	labile_carbon_d = 0.0;
	
	pH = 6.5;
	
	dperc = 0.0;

	solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
	solvesomcent_endyr   = (int)(SOLVESOMCENT_SPINEND   * (nyear_spinup - freenyears) + freenyears);

	temp25 = 10;

	/////////////////////////////////////////////////////
	// Arctic and wetland initialisation

	// Daily heterotrophic respiration. Only ever nonzero on peatland stands.
	dcflux_soil = 0.0;

	// Indices
	IDX = 0;
	SIDX = 0;
	MIDX = 0;
	SIDX_old = 0;

	firstTempCalc = true;
	firstHydrologyCalc = true;
	SIDX_old = 0;
	IDX = 0;
	nsublayer1 = 0;
	nsublayer2 = 0;
	num_evaplayers = 0;

	// initialise snow variables
	snow_active = false;
	snow_active_layers = 0;
	snow_days = 0;
	snow_days_prev = 365;
	dec_snowdepth = 0.0;

	// Peatland hydrology variables
	awtp = 0.0;
	acro_depth = 0.0;
	cato_depth = 0.0;

	acro_co2 = 934.0; // [mimol L-1]

	wtd = 0.0; // [-100, +300] mm

	acro_por = acrotelm_por - Fgas;
	cato_por = catotelm_por - Fgas;

	Wtot = 0.0;
	stand_water = 0.0;

	dmoss_wtp_limit = 1.0;
	dgraminoid_wtp_limit = 1.0;

	k_O2 = 0.0;
	k_CO2 = 0.0;
	k_CH4 = 0.0;
	Ceq_O2 = 0.0;
	Ceq_CO2 = 0.0;
	Ceq_CH4 = 0.0;

	ch4_store = 0.0;
	co2_store = 0.0;

	maxthawdepththisyear = 0.0;
	thaw = 0.0;
	snowdens = snowdens_start; // kg/m3
	dsnowdepth = 0.0;

	for (int mth = 0; mth<12; mth++) {
		mthaw[mth] = 0.0;
		msnowdepth[mth] = 0.0;
		mwtp[mth] = 0.0;

		for (int sl = 0; sl<SOILTEMPOUT; sl++)
			T_soil_monthly[mth][sl] = -999.0;
	}

	for (int sly = 0; sly<NSUBLAYERS_ACRO; sly++)
		sub_water[sly] = 0.0;

	for (int d = 0; d<365; d++) {

		wtp[d] = 0.0;

		for (int ly = 0; ly<NLAYERS; ly++) {

			Frac_ice[ly] = 0.0;
			T_soil[ly] = 0.0;

			if (d == 0) {
				// First day only:
				Frac_water[ly] = 0.0;
				T[ly] = 0.0;
				T_old[ly] = 0.0;
				Frac_org[ly] = 0.0;
				Frac_peat[ly] = 0.0;
				Frac_min[ly] = 0.0;
				Fpwp_ref[ly] = 0.0;
				Frac_water_belowpwp[ly] = 0.0;
				Frac_ice_yesterday[ly] = 0.0;
				Dz[ly] = 0.0;
				por[ly] = 0.0;
				rootfrac[ly] = 0.0;
				CH4_yesterday[ly] = 0.0;
				CO2_soil_yesterday[ly] = 0.0;
				volume_liquid_water[ly] = 0.0;
				total_volume_water[ly] = 0.0;
				tiller_area[ly] = 0.0001; // The max value
				CH4_ebull_ind[ly] = 0.0;
				CH4_ebull_vol[ly] = 0.0;
				CH4_gas_yesterday[ly] = 0.0;
				CH4_diss_yesterday[ly] = 0.0;
				CH4_gas_vol[ly] = 0.0;
			}

			CH4[ly] = 0.0;
			CO2_soil[ly] = 0.0;
			O2[ly] = 0.0;
			CO2_soil_prod[ly] = 0.0;
			CH4_prod[ly] = 0.0;
			CH4_gas[ly] = 0.0;
			CH4_diss[ly] = 0.0;
			Frac_air[ly] = 0.0;
		}
	}

	for (int sl = 0; sl<NSOILLAYER; sl++) {
		wcont[sl] = 0.0;
		alwhc[sl] = 0.0;
		alwhc_init[sl] = 0.0;
	}

	awcont_upper = 0.0;

	nsublayer1 = (int)(SOILDEPTH_UPPER / Dz_soil);		// 5, typically
	nsublayer2 = (int)(SOILDEPTH_LOWER / Dz_soil);		// 10, typically
	num_evaplayers = (int)(SOILDEPTH_EVAP / Dz_soil);	// 2, typically

	// Depths of the acrotelm & catotelm 
	acro_depth = NACROTELM * Dz_acro;
	cato_depth = SOILDEPTH_UPPER + SOILDEPTH_LOWER - acro_depth;

	// end of initialisation
}

/// Get soil clay fraction
double Soil::get_clayfrac() {

	// Return soil clay fractionh
	if (patch.stand.is_highlatitude_peatland_stand())
		return soiltype.clay_frac_peat;
	else
		return soiltype.clay_frac;
}

/// Get soil sand fraction
double Soil::get_sandfrac() {

	// Return soil clay fraction
	if (patch.stand.is_highlatitude_peatland_stand())
		return soiltype.sand_frac_peat;
	else
		return soiltype.sand_frac;
}

/// Get soil silt fraction
double Soil::get_siltfrac() {

	// Return soil silt fraction
	if (patch.stand.is_highlatitude_peatland_stand())
		return soiltype.silt_frac_peat;
	else
		return soiltype.silt_frac;
}

// Get the Water Filled Pore Space of the layer. Currently when Century is restricted to the upper of the old 2 layered LPJ soil column, layer defaults to 0.
double Soil::wfps(int layer = 0) const {
	double wcont = layer == 0 ? get_soil_water_upper() : get_soil_water_lower();

    return (wcont * soiltype.gawc[layer] + soiltype.gwp[layer]) / soiltype.gwsats[layer];
}

// True if phase changes allowed, false otherwise
bool Soil::can_freeze() const {
	return ifallowphasechanges;
}

double Soil::get_soil_temp_25() const {

	// Return soil temperature at 25cm depth

	if (iftwolayersoil)
		return temp_analyticsoln;
	else
		return temp25;
}


/// Performs daily updates to soil temperatures using an analytic function
/**
 * Call each simulation day (if iftwolayersoil true) following update of daily air temperature 
 * prior to canopy exchange and SOM dynamics
 * \param climate   The climate to use to update soil water
 * \param depth     The depth at which the soil temperature is calculated
 */
void Soil::soil_temp_analytic(const Climate& climate, double depth) {

	// DESCRIPTION
	// Calculation of soil temperature at depth (m) (usually middle of upper soil layer).
	// Soil temperatures are assumed to follow surface temperatures according to an
	// annual sinusoidal cycle with damped oscillation about a common mean, and a
	// temporal lag.

	// For a sinusoidal cycle, soil temperature at depth z and time t from beginning
	// of cycle given by (Carslaw & Jaeger 1959; Eqn 52; Jury et al 1991):
	//
	//   (1) t(z,t) = t_av + a*exp(-z/d)*sin(omega*t - z/d)
	//
	//   where
	//     t_av      = average (base) air/soil temperature
	//     a         = amplitude of air temp fluctuation
	//     exp(-z/d) = fractional amplitude of temp fluctuation at soil depth z,
	//                 relative to surface temp fluctuation
	//     z/d       = oscillation lag in angular units at soil depth z
	//     z         = soil depth
	//     d         = sqrt(2*k/omega), damping depth
	//     k         = soil thermal diffusivity
	//     omega     = 2*PI/tau, angular frequency of oscillation (radians)
	//     tau       = oscillation period (365 days)
	//
	// Here we assume a sinusoidal cycle, but estimate soil temperatures based on
	// a lag (z/d, converted from angular units to days) relative to air temperature
	// and damping ( exp(-z/d) ) of soil temperature amplitude relative to air
	// temperature amplitude. A linear model for change in air temperature with time
	// during the last 31 days is used to estimate 'lagged' air temperature. Soil
	// temperature today is thus given by:
	//
	//   (2) temp_soil = atemp_mean + exp( -alag ) * ( temp_lag - temp_mean )
	//
	//   where
	//     atemp_mean = mean of monthly mean temperatures for the last year (deg C)
	//     alag       = oscillation lag in angular units at depth 0.25 m
	//                  (corresponds to z/d in Eqn 1)
	//     temp_lag   = air temperature 'lag' days ago (estimated from linear model)
	//                  where 'lag' = 'alag' converted from angular units to days
	//
	// Soil thermal diffusivity (k) is sensitive to soil water content and is estimated
	// monthly based on mean daily soil water content for the past month, interpolating
	// between estimates for 0, 15% and 100% AWHC (Van Duin 1963; Jury et al 1991,
	// Fig 5.11).


	// conversion factor for soil thermal diffusivity from mm2/s to m2/day
	const double DIFFUS_CONV = 0.0864;
		
	const double HALF_OMEGA = 8.607E-3; // corresponds to omega/2 = pi/365 (Eqn 1)
	// soil depth at which to estimate temperature (m)
	// const double DEPTH = SOILDEPTH_UPPER * 0.0005;
	// conversion factor for oscillation lag from angular units to days (=365/(2*PI))		
	const double LAG_CONV = 58.09;
		
	double a, b; // regression parameters
	double k; // soil thermal diffusivity (m2/day)
	double temp_lag; // air temperature 'lag' days ago (see above; deg C)
	double day[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

	if ((date.year == 0 || date.year == patch.stand.first_year) && date.month == 0 && !date.islastday) {

		// First month of simulation, use air temperature for soil temperature
		temp_analyticsoln = climate.temp;
	}
	else {

		if (date.islastday) {

			// Linearly interpolate soil thermal diffusivity given mean
			// soil water content

			if (mwcontupper < 0.15)
				k = ((soiltype.thermdiff_15 - soiltype.thermdiff_0) / 0.15 *
					mwcontupper + soiltype.thermdiff_0) * DIFFUS_CONV;
			else
				k = ((soiltype.thermdiff_100 - soiltype.thermdiff_15) / 0.85 *
					(mwcontupper - 0.15) + soiltype.thermdiff_15) * DIFFUS_CONV;

			// Calculate parameters alag and exp(-alag) from Eqn 2

			alag = depth / sqrt(k / HALF_OMEGA); // from Eqn 1
			exp_alag = exp(-alag);

		}

		// Every day, calculate linear model for trend in daily air
		// temperatures for the last 31 days: temp_day = a + b * day

		double buffer[31];
		climate.dtemp_31.to_array(buffer);
		regress(day, buffer, 31, a, b);

		// Calculate soil temperature

		temp_lag = a + b * (30.0 - alag * LAG_CONV);

		// Calculate temperatures the old way
		temp_analyticsoln = climate.atemp_mean + exp_alag*(temp_lag - climate.atemp_mean);
	}
}


/// Updates the water content for each soil layer
/**
 * Call each simulation day (if iftwolayersoil true). 
 * Soil class method used to update the water content for each soil layer
 * given snow melt, rainfall, evapotranspiration from vegetation (AET) and
 * percolation between layers. Runoff is also calculated.
 *
 * \param climate   The climate to use to update soil water
 * \param fevap     fraction of modelled area (grid cell or patch) subject to
 *                  evaporation from soil surface
 */
void Soil::hydrology_lpjf_twolayer(const Climate& climate, double fevap) {

	// DESCRIPTION
	// This version replicates the hydrology_lpjf function in LPJ-GUESS v4.0, so 
	// it uses 2 soil layers, and does not consider impedance or ice formation. 

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// UPDATED PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer, 14=bottom layer] as
	//              fraction of available water holding capacity (AWC)
	//				This method ensures that layers 2-4 have identical wcont values, as do layers 5-14
	// wcont_evap = water content of evaporation sublayers (20 cm) at top of soil column
	//              as fraction of available water holding capacity (AWC)
	//				This method ensures that layers 0-1 have identical wcont_evap values
	// awcont_upper = wcont averaged over the growing season
	// dperc      = daily percolation beyond system (mm)
	// runoff     = total daily runoff from all soil layers (mm/day)

	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET) and percolation between layers;
	// calculation of runoff

	// INPUT PARAMETERS
	// rain_melt  = inward water flux to soil today (rain + snowmelt) (mm)
	// perc_base  = coefficient in percolation calculation (K in Eqn 31, Haxeltine
	//              & Prentice 1996)
	// perc_exp   = exponent in percolation calculation (=4 in Eqn 31, Haxeltine &
	//              Prentice 1996)
	// awc        = array containing available water holding capacity of soil
	//              layers (mm rainfall) [0=upper layer]
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface
	// snowpack   = depth of snow (mm)

	// INPUT AND OUTPUT PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer] as
	//              fraction of available water holding capacity (AWC)
	// wcont_evap = water content of evaporation sublayer at top of upper soil layer
	//              as fraction of available water holding capacity (AWC)
	// awcont_upper = wcont averaged over the growing season - guess2008
	// dperc      = daily percolation beyond system (mm)

	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)

	const double SOILDEPTH_EVAP = 200.0;
	// depth of sublayer at top of upper soil layer, from which evaporation is
	// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double BASEFLOW_FRAC = 0.5;
	// Fraction of standard percolation amount from lower soil layer that is
	// diverted to baseflow runoff
	const double K_DEPTH = 0.4;
	const double K_AET = 0.52;
	// Fraction of total (vegetation) AET from upper soil layer that is derived
	// from the top K_DEPTH (fraction) of the upper soil layer
	// (parameters for calculating K_AET_DEPTH below)
	const double K_AET_DEPTH = (SOILDEPTH_UPPER / SOILDEPTH_EVAP - 1.0) *
		(K_AET / K_DEPTH - 1.0) / (1.0 / K_DEPTH - 1.0) + 1.0;
	// Weighting coefficient for AET flux from evaporation layer, assuming active
	//   root density decreases with soil depth
	// Equates to 1.3 given SOILDEPTH_EVAP=200 mm, SOILDEPTH_UPPER=500 mm,
	//   K_DEPTH=0.4, K_AET=0.52

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_water_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	double aet;								// AET for a particular layer and individual (mm)
	double aet_twolayer[2] = {0.0, 0.0};	// total AET for each Gerten-type soil layer (mm)
	double perc_frac;
	double aet_total = 0.0;

	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int s = 0; s<NSOILLAYER; s++) {
			aet = patch.pft[indiv.pft.id].fwuptake[s] * indiv.aet;
			aet_total += aet;
			if (s < NSOILLAYER_UPPER)
				aet_twolayer[0] += aet;
			else
				aet_twolayer[1] += aet;
		}
		vegetation.nextobj();
	}

	// Get the two-layer versions of the wcont values
	double wcont_twolayer[2] = { 0.0, 0.0 };	// wcont for each Gerten-type soil layer
	wcont_twolayer[0] = get_soil_water_upper();
	wcont_twolayer[1] = get_soil_water_lower();
	double twolayerwcont_evap = get_layer_soil_water_evap();

	// Evaporation from soil surface

	// guess2008 - changed to wcont_evap**2, as in LPJ-mL
	// - see Bondeau et al. (2007),  Rost et al. (2008)
	// Added the snowdepth restriction too.
	double evap = 0.0;
	if (snowpack < 10.0) {					// evap only if snow depth < 10mm
		evap = climate.eet * PRIESTLEY_TAYLOR * twolayerwcont_evap * twolayerwcont_evap * fevap;
	}

	// Implement in- and outgoing fluxes to upper soil layer
	// BLARP: water content can become negative, though apparently only very slightly
	//    - quick fix implemented here, should be done better later

	wcont_twolayer[0] += (rain_melt - aet_twolayer[0] - evap) / soiltype.gawc[0];
	if (wcont_twolayer[0] != 0.0 && wcont_twolayer[0] < 0.0001) { // guess2008 - bugfix
		wcont_twolayer[0] = 0.0;
	}

	// Surface runoff
	double runoff_surf = 0.0;
	if (wcont_twolayer[0] > 1.0) {
		runoff_surf = (wcont_twolayer[0] - 1.0) * soiltype.gawc[0];
		wcont_twolayer[0] = 1.0;
	}

	// Update water content in evaporation layer for tomorrow

	twolayerwcont_evap += (rain_melt - aet_twolayer[0] * SOILDEPTH_EVAP*K_AET_DEPTH / SOILDEPTH_UPPER - evap)
		/ soiltype.gawc[0];

	if (twolayerwcont_evap < 0) {
		twolayerwcont_evap = 0;
	}

	if (twolayerwcont_evap > wcont_twolayer[0]) {
		twolayerwcont_evap = wcont_twolayer[0];
	}

	// Percolation from evaporation layer
	double perc = 0.0;
	if (percolate) {
		perc = min(SOILDEPTH_EVAP / SOILDEPTH_UPPER*soiltype.perc_base*pow(twolayerwcont_evap, soiltype.perc_exp),
			max_rain_melt);
	}
	twolayerwcont_evap -= perc / soiltype.gawc[0];

	// Percolation and fluxes to and from lower soil layer(s)

	// Transfer percolation between soil layers
	// Excess water transferred to runoff
	// Eqns 26, 27, 31, Haxeltine & Prentice 1996

	double runoff_drain = 0.0;

	for (int s = 1; s<NSOILLAYER_SIMPLE; s++) {

		// Percolation
		// Allow only on days with rain or snowmelt (Dieter Gerten, 021216)

		if (percolate) {
			perc = min(soiltype.perc_base*pow(wcont_twolayer[s - 1], soiltype.perc_exp), max_rain_melt);
		}
		else {
			perc = 0.0;
		}
		perc_frac = min(perc / soiltype.gawc[s - 1], wcont_twolayer[s - 1]);

		wcont_twolayer[s - 1] -= perc_frac;
		wcont_twolayer[s] += perc_frac * soiltype.gawc[s - 1] / soiltype.gawc[s];
		if (wcont_twolayer[s] > 1.0) {
			runoff_drain += (wcont_twolayer[s] - 1.0)*soiltype.gawc[s];
			wcont_twolayer[s] = 1.0;
		}

		// Deduct AET from this soil layer
		// BLARP! Quick fix here to prevent negative soil water

		wcont_twolayer[s] -= aet_twolayer[s] / soiltype.gawc[s];
		if (wcont_twolayer[s] < 0.0) {
			wcont_twolayer[s] = 0.0;
		}
	}

	// Baseflow runoff (Dieter Gerten 021216) (rain or snowmelt days only)
	double runoff_baseflow = 0.0;
	if (percolate) {
		double perc_baseflow = BASEFLOW_FRAC*soiltype.perc_base*pow(wcont_twolayer[NSOILLAYER_SIMPLE - 1], soiltype.perc_exp);
		// guess2008 - Added "&& rain_melt >= runoff_surf" to guarantee nonnegative baseflow.
		if (perc_baseflow > rain_melt - runoff_surf && rain_melt >= runoff_surf) {
			perc_baseflow = rain_melt - runoff_surf;
		}

		// Deduct from water content of bottom soil layer

		perc_frac = min(perc_baseflow / soiltype.gawc[NSOILLAYER_SIMPLE - 1], wcont_twolayer[NSOILLAYER_SIMPLE - 1]);
		wcont_twolayer[NSOILLAYER_SIMPLE - 1] -= perc_frac;
		runoff_baseflow = perc_frac * soiltype.gawc[NSOILLAYER_SIMPLE - 1];
	}

	// Total runoff
	runoff = runoff_surf + runoff_drain + runoff_baseflow;

	// save percolation from system (needed in leaching())
	dperc = runoff_baseflow + runoff_drain;

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain;
	patch.abaserunoff += runoff_baseflow;
	patch.arunoff += runoff;
	patch.aaet += aet_total;
	patch.aevap += evap;

	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;


	// update awcont_upper
	// Original algorithm by Thomas Hickler

	// Reset awcont_upper on the first day of every year
	if (date.day == 0) {
		awcont_upper = 0.0;
		patch.growingseasondays = 0;
	}
	
	// If it's warm enough for growth, update awcont_upper with this day's wcont
	if (climate.temp > 5.0) {
		awcont_upper += wcont_twolayer[0];
		patch.growingseasondays++;
	}

	// Average awcont_upper on the last day of every year
	if (date.islastday && date.islastmonth) {
		if (patch.growingseasondays > 1)
			awcont_upper /= (double)patch.growingseasondays;
		else
			awcont_upper = 0.0; // This will lead to no establishment for this PFT/species if ifdroughtlimitedestab is TRUE
	}


	// Use the two-layer values calculated above to update the soil properties needed in all NSOILLAYER layers

	// First update wcont for the evaporation layers
	for (int s = 0; s<2; s++) {
		set_layer_soil_water(s, twolayerwcont_evap);
	}

	// Now update wcont for the remaining NSOILLAYER_UPPER layers
	// NB: we want to set the wcont[] value in layers 2-4 to be consistent with twolayerwcont_evap and
	//	wcont_twolayer[0] calculated above. Let wcont_new be this value. Thus:
	//	(2 * twolayerwcont_evap + 3 * wcont_new) / 5 = wcont_twolayer[0], 
	//	giving:
	//	wcont_new = (5 * wcont_twolayer[0] - 2 * twolayerwcont_evap) / 3
	double wcont_new = (5.0 * wcont_twolayer[0] - 2.0 * twolayerwcont_evap) / 3.0;
	wcont_new = min(max(wcont_new, 0.0), 1.0); // Limit to valid values before updating wcont[]

	for (int s = 2; s<NSOILLAYER_UPPER; s++) {
		set_layer_soil_water(s, wcont_new);
	}

	// Next update wcont for the deepest NSOILLAYER_LOWER-NSOILLAYER_UPPER layers
	for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {
		set_layer_soil_water(s, wcont_twolayer[1]);
	}

	// Finally, update wcont_evap, Frac_water etc. in this patch.soil object
	update_soil_water();

}


bool Soil::ice_in_top_layer() {

	const double max_ice_fraction = 0.05; // 5% avoids irrigation restrictions when there are tiny amounts of ice in the soil
	for (int ly = 0; ly < NSOILLAYER_UPPER; ly++) {
		if (Frac_ice[ly + IDX] > max_ice_fraction) return true;
	} 

	return false;
}


/// Updates the water content for each soil layer in non-peatland soils
/**
 * Call each simulation day (if iftwolayersoil false).
 * Soil class method used to update the water content for each soil layer
 * given snow melt, rainfall, evapotranspiration from vegetation (AET) and
 * percolation between layers. Runoff is also calculated.
 *
 * \param climate   The climate to use to update soil water
 * \param fevap     fraction of modelled area (grid cell or patch) subject to
 *                  evaporation from soil surface
 */
void Soil::hydrology_lpjf(const Climate& climate, double fevap) {

	// DESCRIPTION
	// Soil class method used to update the water content for each soil layer 
	// given snow melt, rainfall, evapotranspiration from vegetation (AET) and 
	// percolation between layers. Runoff is also calculated.
	// Updated to deal with 15 soil layers and ice formation.  

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// UPDATED PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer, 14=bottom layer] as
	//              fraction of available water holding capacity (AWC)
	// wcont_evap = water content of evaporation sublayers (20 cm) at top of soil column
	//              as fraction of available water holding capacity (AWC)
	// awcont     = wcont averaged over the growing season
	// dperc      = daily percolation beyond system (mm)
	// runoff     = total daily runoff from all soil layers (mm/day)

	// depth of sublayer at top of upper soil layer, from which evaporation is
	// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double SOILDEPTH_EVAP = 200.0;
	// Fraction of standard percolation amount from lower soil layer that is
	// diverted to baseflow runoff
	const double BASEFLOW_FRAC = 0.5;

	const double MIN_WATER_PERC = 0.0000001;	// mm - Minimum water amount for percolation to occur
	const double maxerr = 0.0001;			// mm - max error allowed

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}


	// *** EVAP - based on initial wcont_evap, as in v4.0 *** 

	double evap_init = 0;
	double Faw_evap_init = 0.0;		// mm
	double awc_init = 0.0;

	for (int s = 0; s < num_evaplayers; s++) {
		Faw_evap_init += wcont[s] * soiltype.awc[s]; // mm;
		awc_init += soiltype.awc[s];
	}

	double wcont_evap_init = Faw_evap_init / awc_init;

	if (snowpack < 10.0) {	// evap only if snow depth < 10mm
							// The potential evaporation:
		evap_init = climate.eet * PRIESTLEY_TAYLOR * wcont_evap_init * wcont_evap_init * fevap;
		// Below, we will limit evaporation to a value <= the available water in the evaporation layer
	}


	// *** AET ***

	// Reset
	double aet = 0.0; // AET for a particular layer and individual (mm)
	double aet_total = 0.0;
	double aet_evap_layers = 0.0;
	double aet_top_layer = 0.0;
	// total AET for each soil layer (mm)
	double aet_layer[NSOILLAYER];
	for (int ly = 0; ly < NSOILLAYER; ly++) {
		aet_layer[ly] = 0.0; // Initialise
	}


	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int ly = 0; ly<NSOILLAYER; ly++) {
			aet = patch.pft[indiv.pft.id].fwuptake[ly] * indiv.aet;
			aet_layer[ly] += aet;
			aet_total += aet;

			// AET from evaporation layers
			if (ly < num_evaplayers)
				aet_evap_layers += aet;

			// AET from top (50cm) layer
			if (ly < NSOILLAYER_UPPER)
				aet_top_layer += aet;
		}
		vegetation.nextobj();
	}


	// *** INITIALISE ***

	// First determine how much water and ice is currently held in each soil layer

	// available water for each soil layer (mm)
	double Faw_layer[NSOILLAYER];
	// ice each soil layer (mm)
	double ice_layer[NSOILLAYER];
	// water that can still be added to each soil layer (mm)
	double potential_layer[NSOILLAYER];
	// Prepare for water balance tests
	double initial_water_in_column = 0.0;


	// Top and bottom layers, initial Faw and potential
	double Faw_layer1 = 0.0;
	double Faw_layer2 = 0.0;
	double potential_layer1 = 0.0;
	double potential_layer2 = 0.0;


	for (int ly = 0; ly < NSOILLAYER; ly++) {

		if (wcont[ly] < 0.0 || wcont[ly] > 1.0) {
			dprintf("Soil::hydrology_lpjf - bad wcont!\n");
			return;
		}

		Faw_layer[ly] = wcont[ly] * soiltype.awc[ly]; // mm

		ice_layer[ly] = Frac_ice[ly + IDX] * Dz[ly + IDX]; // mm

		// Inital water in this layer, before AET
		double layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

		// for water balance checks
		initial_water_in_column += layerwater;

		// Shouldn't happen
		if (aet_layer[ly] > Faw_layer[ly]) {
			int test = 1;
		}

		// AW
		// Now remove AET as it's already accounted for
		Faw_layer[ly] -= min(aet_layer[ly], Faw_layer[ly]);
		// Update layer water after subtracting AET
		layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

		// Update wcont
		wcont[ly] = Faw_layer[ly] / soiltype.awc[ly]; // mm

		if (ly < NSOILLAYER_UPPER)
			Faw_layer1 += Faw_layer[ly]; // mm
		else
			Faw_layer2 += Faw_layer[ly]; // mm


		// POTENTIAL
		potential_layer[ly] = aw_max[ly] - layerwater;

		if (ly < NSOILLAYER_UPPER)
			potential_layer1 += potential_layer[ly]; // mm
		else
			potential_layer2 += potential_layer[ly]; // mm


		// Check balance
		if (potential_layer[ly] < -1.0 * maxerr) {
			// fail("Soil::hydrology_lpjf - error in a soil layer's water balance!\n");	
			dprintf("ERROR: Soil::hydrology_lpjf - error in a soil layer's water balance! potential_layer#: %g\n", ly);
		}
		else {
			potential_layer[ly] = max(0.0, potential_layer[ly]); // catches small negative values
		}
	} // for loop (ly) 


	// *** EVAPORATION and TOP LAYER ***

	// Update wcont_evap
	// wcont_evap is the available, nonfrozen water down to SOILDEPTH_EVAP, divided by the capacity of those layers to hold water 
	wcont_evap = 0.0;
	double awc_count = 0.0;
	double Faw_layer_evap = 0.0;		// mm
	double potential_evap_layer = 0.0;	// mm
	double potential_top_layer = 0.0;	// mm

	for (int s = 0; s < num_evaplayers; s++) {
		Faw_layer_evap += Faw_layer[s];
		awc_count += soiltype.awc[s];
		potential_evap_layer += potential_layer[s];
	}

	wcont_evap = Faw_layer_evap / awc_count;


	// *** EVAPORATION FROM BARE SURFACE (after AET removed) ***
	double evap = min(evap_init, max(0.0, Faw_layer_evap));

	// REMOVAL of water from the evaporation layers of the soil (evap > 0).
	if (evap > 0.0 && Faw_layer_evap > 0.0) {

		// Specifying evap > 0 ensures that there is at least some water in the evap layers before trying to remove water

		for (int s = 0; s<num_evaplayers; s++) {

			// Remove the water in proportional to the available water remaining in each layer
			double proportion = evap * (Faw_layer[s] / Faw_layer_evap); // mm
			Faw_layer[s] -= proportion;

			potential_layer[s] += proportion; // potential INCREASES
			wcont[s] = Faw_layer[s] / soiltype.awc[s];
			oob_check_wcont(wcont[s]);
		}
	} // > 0.0?
	else {
		evap = 0.0; // No evaporation in this case
	}


	// Implement in- and outgoing fluxes to the top 50cm soil layers, taking into account AET in those layers and evaporation
	double water_flux_in = rain_melt;	// mm

	// Initialise runoff variables

	// Surface runoff (mm) 
	double runoff_surf = 0.0;
	// Baseflow runoff (Dieter Gerten 021216) (rain or snowmelt days only)
	double runoff_drain = 0.0;
	double runoff_baseflow = 0.0;


	// Potential in the top layer
	potential_top_layer = 0.0;
	for (int s = 0; s < NSOILLAYER_UPPER; s++) {
		potential_top_layer += potential_layer[s];
	}

	// *** INPUT TO TOP LAYER ***

	if (water_flux_in > 0) {

		// ADDITION of water to the top layers (rain_melt > 0)
		if (water_flux_in < potential_top_layer) {

			// The upper soil layers can absorb all of this water today
			// so add the water to layers in proportion to their capacity
			for (int s = 0; s<NSOILLAYER_UPPER; s++) {

				double water_input_ly = 0.0;
				if (potential_top_layer > 0.0)
					water_input_ly = water_flux_in * (potential_layer[s] / potential_top_layer);

				Faw_layer[s] += water_input_ly;
				potential_layer[s] -= water_input_ly;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);
			}

			// No surface runoff
			runoff_surf = 0.0;
		}
		else {

			// The upper soil layers cannot absorb all of this water today, so take up what can be absorbed 
			// (i.e. potential_top_layer) and add the rest to surface runoff

			for (int s = 0; s<NSOILLAYER_UPPER; s++) {
				Faw_layer[s] += potential_layer[s];
				potential_layer[s] = 0.0;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
			}

			// Update surface runoff
			runoff_surf = water_flux_in - potential_top_layer;
		}

		for (int s = 0; s<NSOILLAYER_UPPER; s++) {
			oob_check_wcont(wcont[s]);
		}
	}


	// Update wcont_evap again, after the input of rain_melt
	wcont_evap = 0.0;
	Faw_layer_evap = 0.0;	// mm

	for (int s = 0; s < num_evaplayers; s++) {
		Faw_layer_evap += Faw_layer[s];
	}

	wcont_evap = Faw_layer_evap / awc_count;
	oob_check_wcont(wcont_evap);


	// Update Faw_layer1 and wcont1

	double awc1 = 0.0;
	Faw_layer1 = 0.0;
	for (int s = 0; s < NSOILLAYER_UPPER; s++) {
		Faw_layer1 += Faw_layer[s];
		awc1 += soiltype.awc[s];
	}

	double wcont_layer1 = Faw_layer1 / awc1;


	// *** PERCOLATION ***
	// We assume there are two percolation layers. 
	// 1: the top 50cm soil layers considered in the Gerten et al. scheme
	// 2: Lower 1.0 m layer considered in the Gerten et al. scheme

	if (percolate) {

		// Percolation from top layer if there is enough water (at least MIN_WATER_PERC)
		double perc = 0.0;

		if (Faw_layer1 > MIN_WATER_PERC) {

			// Percolation from the top  layer (limited to available liquid water)
			perc = min(soiltype.perc_base * pow(wcont_layer1, soiltype.perc_exp), Faw_layer1);

			for (int s = 0; s<NSOILLAYER_UPPER; s++) {

				// Take water proportionally from the layers.
				double perc_layer = perc * (Faw_layer[s] / Faw_layer1);
				Faw_layer[s] -= perc_layer;
				potential_layer[s] += perc_layer;

				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);
			}
		}


		// This is the percolation (mm) from the top layers
		double perc_from_above = perc;

		// Percolation and fluxes to and from the soil layers below the top 50cm

		// Excess water transferred to runoff
		// Eqns 26, 27, 31, Haxeltine & Prentice 1996

		double water_input = perc_from_above;

		// First we need to add the water percolated from above (perc_from_above) to this layer
		if (perc_from_above > potential_layer2) {
			water_input = potential_layer2; // Fill the layer and reduce water_input
			// Update runoff_drain
			runoff_drain = perc_from_above - potential_layer2; // The excess is run off and not percolated further
		}

		perc_from_above = 0.0;

		double Faw_perc_layer_2 = 0.0;
		double awc_perc_layer_2 = 0.0;
		double wcont_perc_layer_2 = 0.0;

		if (water_input >= 0) {

			// Now add any water (water_input) to layers 5-14 in proportion to the capacity, and recalculate Faw, wcont and awc
			for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {

				double water_input_ly = 0.0;
				if (potential_layer2 > 0.0)
					water_input_ly = water_input * (potential_layer[s] / potential_layer2);

				Faw_layer[s] += water_input_ly;
				potential_layer[s] -= water_input_ly;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);

 				Faw_perc_layer_2 += Faw_layer[s];
				awc_perc_layer_2 += soiltype.awc[s];
			}

			// Update wcont for these layers before percolation
			wcont_perc_layer_2 = Faw_perc_layer_2 / awc_perc_layer_2;
		}

		// Percolation 3 -> base

		double perc_from_base = 0.0;
		if (percolate && Faw_perc_layer_2 > MIN_WATER_PERC) {

			perc_from_base = min(BASEFLOW_FRAC * soiltype.perc_base * pow(wcont_perc_layer_2, soiltype.perc_exp), max_rain_melt);
			perc_from_base = min(perc_from_base, Faw_perc_layer_2); // Only available water in these layers can percolate

			// As in LPJ-GUESS v4.0
			if (perc_from_base > rain_melt - runoff_surf && rain_melt >= runoff_surf)
				perc_from_base = rain_melt - runoff_surf;
		}
		else {
			// No percolation
			perc_from_base = 0.0;
		}

		// Now remove the water (i.e. perc_from_base) from layers 5-14 in proportion to the available water
		if (perc_from_base > 0) {

			for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {

				// Remove the water (perc_from_base) in proportional to the available water remaining in each layer
				double proportion = perc_from_base * (Faw_layer[s] / Faw_perc_layer_2); // mm

				Faw_layer[s] -= proportion;

				potential_layer[s] += proportion;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);
			}
		}

		// Baseflow
		runoff_baseflow = perc_from_base;
	}
	else {

		// No percolation at all, so reset runoff to 0 
		runoff_baseflow = 0.0;
		runoff_drain = 0.0;
	}


	// Total runoff
	runoff = runoff_surf + runoff_drain + runoff_baseflow;

	// water added when patch.stand.is_true_wetland_stand() should be be subtracted from runoff
	// in proportion to its components
	if (patch.stand.is_true_wetland_stand() && ifsaturatewetlands) {

		if (runoff <= patch.wetland_water_added_today && runoff > 0.0) {

			// Not enough runoff to balance the water added to the wetland, i.e. we added more than was run off.

			// Update wetland_water_added by taking back the amount used for runoff.
			patch.wetland_water_added_today -= runoff;

			// All runoff removed 
			runoff_surf = 0.0;
			runoff_baseflow = 0.0;
			runoff_drain = 0.0;
			runoff = 0.0;

		}
		else if (runoff > patch.wetland_water_added_today  && runoff > 0.0) {

			// Enough runoff to balance the water added to the wetland so we take it back.
			runoff_surf -= patch.wetland_water_added_today * runoff_surf / runoff;
			runoff_drain -= patch.wetland_water_added_today * runoff_drain / runoff;
			runoff_baseflow -= patch.wetland_water_added_today * runoff_baseflow / runoff;
			// Recalculate total runoff 
			runoff = runoff_surf + runoff_drain + runoff_baseflow;

			// Reset wetland_water_added
			patch.wetland_water_added_today = 0.0;

		}
	}

	// save percolation from system (needed in leaching())
	dperc = runoff_baseflow + runoff_drain;

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain;
	patch.abaserunoff += runoff_baseflow;
	patch.arunoff += runoff;
	patch.awetland_water_added += patch.wetland_water_added_today;
	patch.aaet += aet_total;
	patch.aevap += evap;

	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;

	// *** WATER IN BALANCE? ***
	if (DEBUG_SOIL_WATER) {

		double final_water_in_column = 0.0;

		for (int ly = 0; ly < NSOILLAYER; ly++) {

			Faw_layer[ly] = wcont[ly] * soiltype.awc[ly]; // mm
			ice_layer[ly] = Frac_ice[ly + IDX] * Dz[ly + IDX]; // mm

			// Water in this layer
			double layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

			// for water balance checks
			final_water_in_column += layerwater;

		} // for loop (ly)

		  // is water in + initial storage = water out + final storage?
		double water_in_storage_in = initial_water_in_column + rain_melt;
		double water_out_storage_out = final_water_in_column + evap + aet_total + runoff;

		if (fabs(water_in_storage_in - water_out_storage_out) > maxerr) {
			dprintf("Soil::hydrology_lpjf - error in the TOTAL water balance!\n");
			return;
		}
	}

	// Drought limited establishment - update awcont_upper
	// Original algorithm by Thomas Hickler
	for (int s = 0; s<NSOILLAYER_UPPER; s++) {

		// Reset on the first day of every year
		if (date.day == 0 && s == 0) {
			awcont_upper = 0.0;
			patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			awcont_upper += wcont[s] / (double)NSOILLAYER_UPPER;
			if (s == 0) {
				patch.growingseasondays++;
			}
		}
	}

	// Average awcont_upper on the last day of every year
	if (date.islastday && date.islastmonth) {
		if (patch.growingseasondays > 1)
			awcont_upper /= (double)patch.growingseasondays;
		else
			awcont_upper = 0.0; // This will lead to no establishment for this PFT/species if ifdroughtlimitedestab is TRUE
	}

	// Finally, recalculate Frac_water based on these calculations and updates to wcont[] 
	update_soil_water();
}


/// Updates the water content for each soil layer in peatland soils
/**
* Call each simulation day (if iftwolayersoil false).
* Soil class method used to update the water content for each soil layer
* given snow melt, rainfall, evapotranspiration from vegetation (AET) and
* percolation between layers. Runoff is also calculated.
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soil::hydrology_peat(const Climate& climate, double fevap) {

	// DESCRIPTION
	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET) and percolation between layers;
	// calculation of runoff
	// Uses peatland hydrology scheme from Granberg et al. (1999) and Wania et al. (2009a, 2009b)

	// INPUT PARAMETERS
	// climate    = climate today (mm)
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface


	double aet_total = 0.0;		// total AET (mm)
	double runoff_surf;			// runoff from upper soil layer (mm)
	double runoff_drain;		// runoff (drainage) from lower soil layers (mm)
	double evap;				// evaporation from soil surface (mm)
	double evapotranspiration;	// evapotranspiration (mm/day)


	// PARAMETERS and variables for the calculation of the WTP
    // See: Granberg et al. (1999), Eqns. 1-5.
	
	/// depth of each sphagnum layer [mm]
	const double Dz_sph = 10.0;
	/// surface water content when WTD is at 20 mm
	const double surfw = 0.34;
	/// minimum fractional water content at surface in mm3/mm3 - see Wania et al (2009a), Eqn 22. Granberg et al. (1999), Eqns. 1-3.
	const double minvtot = 0.25;              
	/// The soil water characteristics in the upper peat are considered to be linear in the suction interval 0-zmin and constant 
	/// in the interval zmin - acrotelm depth (normally 300mm). See: Granberg et al. (1999), Eqns. 1-3.
	const double zmin = 100.0;
	/// gradient limit of upper gradient [mm]
	const double grad_lim = 50.0;
	/// Granberg et al. (1999) parameter, eqn 2.
	double az; 
	/// Depth of sublayers in the acrotelm
	double Dz_sub; 	

	/// Day of year
	int daynum = date.day;

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_water_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	// Update Wtot, whc[], awhc[], and wcont[]
	if (!update_layer_water_content(daynum)) 
		fail();

	// Porosity updates in the presence of ice? 
	bool subtractIceFromWtot = false;
	double acro_por_icy = acro_por; // the pore space minus the ice fraction
	double avgIceFrac = 0.0;
	
	if (subtractIceFromWtot) {	
		for (int ly = IDX; ly < IDX + NACROTELM; ly++) 
			avgIceFrac += (Frac_ice[ly] + Fpwp_ref[ly] - Frac_water_belowpwp[ly])/(double)NACROTELM;
	}

	acro_por_icy -= avgIceFrac;


	// *** TRANSPIRATION *** 
	
	// Retrieve Vegetation for this patch
	Vegetation& vegetation=patch.vegetation;

	// Transpiration from the acrotelm  
	double aet_acrotelm = 0.0;
	double aet = 0.0; // AET for a particular layer and individual (mm)

	// Sum AET for across all vegetation individuals in this patch

	// Calculate fpc_moss
	double fpc_moss_total = 0.0;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		double pft_acro_root_frac = 0.0;

		for (int ly = IDX; ly < IDX + NACROTELM; ly++)
			pft_acro_root_frac += patch.pft[indiv.pft.id].pft.rootdist[ly-IDX];

		// this includes mosses, which have root_frac = 1 in the acrotelm
		// Fraction of the aet from the acrotelm only.
		aet_acrotelm += indiv.aet * pft_acro_root_frac;
		
		// Total aet, assumed from acrotelm and catotelm combined 
		aet_total += indiv.aet;

		if (indiv.pft.ismoss())
			fpc_moss_total += indiv.fpc;

		vegetation.nextobj();
	}


	// *** EVAPORATION *** 
	// Evaporation from soil surface

	// Wania et al (2009)
	if (snowpack < 10.0 && wtd < SOILDEPTH_EVAP && !negligible(fevap)) { // i.e. evap only if snow depth < 1cm
		evap = fevap*climate.eet*PRIESTLEY_TAYLOR*(0.99 / (1+exp(-1.0*(-wtd + 98.7)/22.6) + 0.02));
	} 
	else {
		evap = 0.0;
	}
	
	// *** EVAPOTRANSPIRATION *** 

	bool empirical_evapotranspiration = false;
	if (empirical_evapotranspiration) {

		// An empirical evapotranspiration parameterisation from Kim & Verma

		/*
		// From LPJ-WHyMe:
		Calculate transpiration and evaporation from wetlands together.
		dpet is the equilibrium ET but is the same as the potential ET in
		Kim and Verma (1996).  (BAD choice of variable names in original
		LPJ!)  Hence, we can use dpet as approximation for potential
		evapotranspiration and scale it by the water table depth.  daept
		is allowed to happen even if the wtp is at -300mm because plants
		will be able to extract water from the unsaturated acrotelm or the
		saturated catotelm.  This however makes it impossible to close the
		hydrological cycle in LPJ-WHy. 
		Kim and Verma (1996): ET = PET * (1.02 + 0.00075 * WTP)
		*/ 

		// Note 1. that this is now the TOTAL evapotranspiration, replacing aet_total

		// Note 2: An alternative parameterisation was given by Yurova et al. (2006)
		// evapotranspiration=eet*exp[0.02*wtp_cm], where // wtp_cm [100,-30]
		// This would give almost identical results to the Kim & Verma formula below 
		// when the water table is at the surface

		if (snowpack < 10.0) { // i.e. evap only if snow depth < 1cm
			evapotranspiration = climate.eet*(1.02 - 0.00075 * wtd)*fevap; // wtd [-100,+300]
		}
		else {
			evapotranspiration = 0.0;
		}

	} 
	else {

		// This is JUST evapotranspiration from the acrotelm
		evapotranspiration = aet_acrotelm + evap;
	}

	// *** ACROTELM WATER VOLUME ***
	
	az = (acro_por - minvtot) / zmin; // Granberg (1999) - Eqn 2 & Wania et al. (2009a), Eqn 22
	
	// Drainage - could possibly be read in for site-specific studies
	runoff_drain = 0; 

	// Update available water 
	Wtot += rain_melt-evapotranspiration-runoff_drain;
	 
	// *** RUNOFF AND RUNON***

	double acrowater = 0.0; // Frac_water[MIDX] * Dz[MIDX];

	for (int ly = IDX; ly < IDX + NACROTELM; ly++)
		acrowater += (Frac_water[ly] + Frac_water_belowpwp[ly]) * Dz[ly];

	double ideal_runoff = exp(-0.01 * wtd);
	if (acrowater > 0.0 && Frac_ice[IDX] < 0.7) // Granberg et al. (1999) use the same ice condition - see f_icestop in their Table 1
		runoff_surf = min(ideal_runoff,acrowater);
	else
		runoff_surf = 0.0;

	// .ins file option: 
	// run on or off - but only when there's flowing water

	double runon = 0.0;

	// runon set in global.ins
	soiltype.runon = wetland_runon;

	if (runoff_surf > 0.0) {

		double runofforon = soiltype.runon;

		if (runofforon < 0) 
			runoff_surf -= runofforon;	// e.g. 3 mm/day for bog-like/hummock conditions

		if (runofforon > 0)
			runon += runofforon;		// 3-5 mm/day for fen-like conditions
	}

	double Wtot_init = Wtot;

	Wtot += runon - runoff_surf;

	// Add standing water, if any
	Wtot += stand_water;

	// Not in LPJ-WHy
	if (Wtot < 0.0) {
		Wtot = 0.0;
		runoff_surf = Wtot_init;
	}

	// Save evapotranspiration data
	if (empirical_evapotranspiration) {

		patch.aaet+=(1.0-fevap)*evapotranspiration;
		patch.aevap+=fevap*evapotranspiration;
		patch.maet[date.month]+=(1.0-fevap)*evapotranspiration;
		patch.mevap[date.month]+=fevap*evapotranspiration;
	} 
	else {
	
		patch.aevap+=evap;
		patch.mevap[date.month]+=evap;
		patch.aaet += aet_total;
		patch.maet[date.month] += aet_total;
	}


	// *** UPDATE WATER TABLE DEPTH ***

	// LPJ-WHy comments:
	// If total water volume is greater than the porosity then the
	// WTP is above the surface. The maximum height of standing water
	// is limited by maxh. If the watertable is above that height, it
	// will run off.  Snow melt season: 1. If there is a standing
	// water layer, which is obviously frozen in winter, then the
	// snow melt will run off and only when all the ice of the stand
	// layer is thawed, the water will be added to the soil
	// layers. 2. If there is no standing water layer, then the
	// snowmelt water...should penetrate into the soil.
	
	// NB 
	// e.g. Wtot = 270mm =>
	// Frac_water = 270/300 = 0.9 in the three acrotelm layers
	
	// Wtot <= 140.03611111 =>
	// Frac_water = 0.27419, 0.4186, 0.7075 in the top three layers (if no ice)

	if (Wtot > acro_depth * acro_por) {
		
		// Saturated acrotelm

		wtp[daynum] = Wtot - acro_depth * acro_por; 
		// i.e. wtp[daynum] > 0 above the surface
		
		if (wtp[daynum] > maxh) {
			runoff_surf += wtp[daynum] - maxh; 
			// Update runoff from the surface
			wtp[daynum] = maxh;
		}

		wtd = -wtp[daynum]; // i.e. wtd defined as < 0 above the surface
		Wtot = acro_depth * acro_por; 
		// I.e. Wtot does NOT include standing water
		stand_water = 0.0; // = wtp[daynum]; // we assume no standing water 
	} 
	else {
		
		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (acro_por*acro_depth - Wtot)/(2.0 * az)); 
		
		if (wtd > zmin)
			wtd = 3.0 * (acro_por * acro_depth - Wtot) / (2.0 * (acro_por - minvtot));

		if (wtd > acro_depth) 
			wtd = acro_depth; 

		wtp[daynum] = -wtd;
		stand_water = 0.0;
	}


	// *** RUNOFF ***

	// Save today's runoff in Soil::runoff
	runoff = runoff_surf;
	patch.arunoff+=runoff;
	patch.mrunoff[date.month]+=runoff;
	patch.asurfrunoff += runoff;


	// Depth of acrotelm sublayers - usually 10 mm
	Dz_sub = NACROTELM * Dz_acro / NSUBLAYERS_ACRO;


	// *** WATER CONTENT IN EACH LAYER - THE PROFILE ***

	double surfw_sat = 0.0;
	double value[NSUBLAYERS_ACRO];
	double acro_aw = 0.0; // mm of plant available water in acrotelm
	double cato_aw = 0.0; // mm of plant available water in catotelm

	if (wtd > 0.0) {

		// Below the surface

		// Determine the subsurface layer (0 to NSUBLAYERS_ACRO) within which wtd lies. 
		// Possible values:
		// = 0 if wtd < 10 mm 
		// = NSUBLAYERS_ACRO iff wtd = Dz_acro
		int wtd_layer = (int)(wtd / Dz_sub); // Takes the integer part only

		surfw_sat = max(minvtot, acro_por - az * wtd); // Granberg - eqn 3.

		if (wtd_layer == 0) {
			
			// Total saturation throughout the acrotelm (from 10mm down)
			for (int j = 0; j < NSUBLAYERS_ACRO; j++) {
				value[j] = acro_por; 
			}
		} 
		else {
			
			// Unsaturated zone in the acrotelm - use Granberg's eqn 1.
			for (int i = 0; i < wtd_layer; i++) {
				double z = (i+1)*Dz_sub; // mm
				value[i] = min(acro_por, surfw_sat + (acro_por_icy - surfw_sat) * pow(z/wtd,2));
			}

			// Saturated zone in the acrotelm
			for (int j = wtd_layer; j < NSUBLAYERS_ACRO; j++) {
				value[j] = acro_por; 
			}
		}

		// Now calculate the volumetric water content throughout the acrotelm
		sub_water[0] = 0.5 * (surfw_sat + value[0]);
		for (int j = 1; j < NSUBLAYERS_ACRO; j++) {
			sub_water[j] = 0.5 * (value[j] + value[j-1]); 
		}

		// Water content in each acrotelm layer
		double Ftotal[NLAYERS]; // Volumetric fraction in the acrotelm layers
		for (int ly = IDX; ly < IDX + NACROTELM; ly++)
			Ftotal[ly] = 0.0;

		int sublayers_per_layer = (int)Dz_acro/(int)Dz_sub; 

		int count = 0;
		for (int ly = IDX; ly < IDX + NACROTELM; ly++) {
		
			double layerwater = 0.0;
			for (int subly = count*sublayers_per_layer; subly < (count+1)*sublayers_per_layer; subly++) {
				layerwater += sub_water[subly] * Dz_sub; // add mm of water for this sublayer 
			}
			
			// Divide by the depth of the acrotelm layers to get the volumetric fraction
			Ftotal[ly] = layerwater / Dz_acro; 
			count++;
		}

		// Subtract ice content, water content below the permanent wilting point
		
		// ACROTELM

		for (int ly = IDX; ly < IDX + NACROTELM; ly++) {
			Frac_water[ly] = max(Ftotal[ly] - Frac_ice[ly] - Fpwp_ref[ly], 0.0);
			acro_aw += Frac_water[ly] * Dz[ly];
		}
	} 
	else { 

		// WTD <= 0.0
		// at or above the surface - standing water

		// Subtract ice content, water content below the permanent wilting point
		// For the calculations below, note that:
		// totaliceinlayer = Frac_ice + Fpwp_ref - Frac_water_belowpwp;
		// totalwaterinlayer = Frac_water + Frac_water_belowpwp;

		// ACROTELM - saturated
		for (int ly = IDX; ly < IDX + NACROTELM; ly++) {
			Frac_water[ly] = max(por[ly]-Frac_ice[ly]-Fpwp_ref[ly], 0.0);
			acro_aw += Frac_water[ly] * Dz[ly];
		}
	}

	// CATOTELM - assumed to be always saturated
	for (int ly = IDX + NACROTELM; ly < NLAYERS; ly++) {
		Frac_water[ly] = max(por[ly] - Frac_ice[ly] - Fpwp_ref[ly], 0.0);
		cato_aw += Frac_water[ly] * Dz[ly];
	}

	// Tidy up Frac_water
	for (int ly = 0; ly < NLAYERS; ly++) {
		if (Frac_water[ly] < Fpwp_ref[ly]) 
			Frac_water[ly] = 0.0;
	}

	// *** WATER CONTENT IN EACH SOIL LATER (FOR PLANT GROWTH) ***
	if (!update_layer_water_content(daynum)) 
		return;

	// UPDATE CO2 LEVEL IN ACROTELM
	update_acrotelm_co2(climate.co2);


	// *** UPDATE LIMITS ON MOSS PHOTOSYNTHESIS ***

	/* 
	The cap is kept at 100% while WTP is above wtp_moss_upper ((Wania et al 2009b) have -150 mm), decreases linearly
	between wtp_moss_upper and -280mm  This has to do with the water
	content of the mosses as they are less able to photosynthesize
	as they dry out. NO MOSS DESSICATION
	*/
	
	double wtp_moss_upper = 0.0;	// Mosses can be too competitive when set to -150 mm, as in Wania et al. 2009b 
	double wtp_moss_lower = -280.0; // Wania et al. 2009b
	double moss_wtp_limit_lowerlimit = 0.3; 

	if (wtp[daynum] > wtp_moss_upper)
		dmoss_wtp_limit = 1.0; // No limit
	else if (wtp[daynum] > wtp_moss_lower)
		dmoss_wtp_limit = 1.0 - (-wtp[daynum] + wtp_moss_upper) * (1.0 - moss_wtp_limit_lowerlimit)/(wtp_moss_upper - wtp_moss_lower);
	else 
		dmoss_wtp_limit = moss_wtp_limit_lowerlimit; // between -300mm and -280mm
	

	// *** UPDATE LIMITS ON GRAMINOID PHOTOSYNTHESIS ***

	// LPJ-WHy code & comment:

	/*
	    IF (pft .EQ. 11 .AND. mwtp(m) .LT. -100) THEN
			! The following lines will reduce the mgpp of
			! flood-tolerant grasses when the water table drops
			! below -10cm. This should reflect that these grasses
			! are usually found where the water table stays close
			! to the surface for most of the time (i.e. in fens).
			! Note that mwtp is negative!
	
			mgpp(m,pft) = mgpp(m,pft) + mgpp(m,pft) *
			$                (mwtp(m) + 100) / 200 
	     ENDIF
	*/

	double wtp_graminoid_upper = -10; // -50.0;		// mm - Wania et al. 2009b have -100
	double wtp_graminoid_lower = -100; // -200.0;	// mm - Wania et al. 2009b have -300 mm
	double graminoid_wtp_limit_lowerlimit = 0.0; 

	// More restrictive than Wania et al., otherwise graminoids become too dominant as they never suffer from inundation stress
	// The cap drops from 1 as soon as wtp goes below wtp_graminoid_upper
	if (wtp[daynum] >= wtp_graminoid_upper)
		dgraminoid_wtp_limit = 1.0; // i.e. no limit!
	else if (wtp[daynum] > wtp_graminoid_lower)
		dgraminoid_wtp_limit = 1.0 - (-wtp[daynum] + wtp_graminoid_upper) * (1.0 - graminoid_wtp_limit_lowerlimit)/(wtp_graminoid_upper - wtp_graminoid_lower);
	else 
		dgraminoid_wtp_limit = graminoid_wtp_limit_lowerlimit; // between -300mm and -200mm


	// *** UPDATE OUTPUT VARIABLES ***

	// Reset monthly wtp and mmoss_wtp_limit averages on Jan 1
	if (daynum == 0) {
		for (int mth = 0; mth < 12; mth++) {	
			mwtp[mth] = 0.0;
		}
	}

	int mdays = date.ndaymonth[date.month];
	mwtp[date.month] += wtp[date.day] / (double)mdays;

	if (firstHydrologyCalc) 
		firstHydrologyCalc = false;

	if (date.day == Date::MAX_YEAR_LENGTH) {

		// Calculate annual average WTP. Needed in update_acrotelm_co2
		awtp = 0.0;
		for (int d = 0; d < 365; d++)
			awtp += wtp[d]/365.0; 
	}

	return; // ...if no problems
}


// return wcont for a certain layer
double Soil::get_layer_soil_water(int layer) const {

	if (layer < 0 || layer > NSOILLAYER-1 || layer != int(layer))
		return -999.0;
	else
		return wcont[layer];
}


// return wcont for a subset of soil layers
double Soil::get_soil_water(int layer1, int layer2) const {

	double total_available_water = 0.0; // available water for the thick soil layer (mm)
	double total_capacity = 0.0;

	if (layer1 < 0 || layer1 > NSOILLAYER || layer1 != int(layer1)) {
		fail("Soil::get_soil_water - bad layer1!\n");
		return -999;
	}

	if (layer2 < 0 || layer2 > NSOILLAYER || layer2 != int(layer2)) {
		fail("Soil::get_soil_water - bad layer2!\n");
		return -999;
	}

	if (layer2 <= layer1) {
		fail("Soil::get_soil_water - layer2 <= layer1\n");
		return -999;
	}

	for (int ly = layer1; ly < layer2; ly++) {
		
		if (wcont[ly] < 0.0 || (wcont[ly] > 1.0 && !negligible(wcont[ly] - 1.0, -12))) {
			fail("Soil::get_soil_water - bad wcont!\n");
			return -999;
		}

		total_available_water += wcont[ly] * soiltype.awc[ly]; // mm
		total_capacity += soiltype.awc[ly]; // mm
	}

	return total_available_water/ total_capacity;
}


// return wcont for the upper 50cm (upper true) or lower (100cm) soil layer (upper false)
double Soil::get_soil_water_upper() const {
	return get_soil_water(0, nsublayer1);
}

// return wcont for the lower (100cm) soil layer
double Soil::get_soil_water_lower() const {

	return get_soil_water(nsublayer1, NSOILLAYER);
}

// return copy of wcont
void Soil::copy_layer_soil_water_array(double wconttoreturn[NSOILLAYER]) {
		
	for (int ly = 0; ly < NSOILLAYER; ly++)
		wconttoreturn[ly] = wcont[ly];
}

// return wcont_evap
double Soil::get_layer_soil_water_evap() const {
		return wcont_evap;
}


// add water to a certain soil layer and return the water that cannot be added due to limited capacity
double Soil::add_layer_soil_water(int layer, double extrawater) {

	double overflow = 0.0;

	if (layer < 0 || layer > NSOILLAYER - 1 || layer != int(layer))
		fail("Soil::add_layer_soil_water - bad layer\n");
	else {

		double Faw_layer = 0.0; // mm
		double ice_layer = 0.0; // mm
		double layerwater_max = 0.0; // mm

		if (patch.stand.is_highlatitude_peatland_stand()) {
			Faw_layer = wcont[layer] * soiltype.awc_peat[layer]; // mm
			ice_layer = Frac_ice[layer + IDX] * Dz[layer + IDX]; // mm
			layerwater_max = (acro_por - peat_wp) * Dz[layer + IDX];

		} 
		else {
			Faw_layer = wcont[layer] * soiltype.awc[layer]; // mm
			ice_layer = Frac_ice[layer + IDX] * Dz[layer + IDX]; // mm
			layerwater_max = aw_max[layer];
		}

		double layerwater = Faw_layer + ice_layer; // mm
		double potential_layer = layerwater_max - layerwater;

		if (extrawater < potential_layer)
			overflow = 0.0;
		else {
			overflow = extrawater - potential_layer;
			extrawater = potential_layer;
		}

		// Now add as much water to this layer as it can hold
		Faw_layer += extrawater;
		potential_layer -= extrawater;

		// wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
		if (patch.stand.is_highlatitude_peatland_stand())
			wcont[layer] = Faw_layer / soiltype.awc_peat[layer];
		else
			wcont[layer] = Faw_layer / soiltype.awc[layer];

		oob_check_wcont(wcont[layer]);
		Frac_water[layer + IDX] = Faw_layer / Dz[layer + IDX];
	}

	return overflow;
}

// set wcont for a certain layer
void Soil::set_layer_soil_water(int layer, double newlayerwater) {

	if (layer < 0 || layer > NSOILLAYER-1 || layer != int(layer))
		fail("Soil::set_layer_soil_water - bad layer\n");
	else {
		wcont[layer] = newlayerwater;
	}
}

// set wcont_evap
void Soil::set_layer_soil_water_evap(double newevaplayerwater) {
		wcont_evap = newevaplayerwater;
}


void Soil::update_soil_water() {

	// DESCRIPTION
	// Update Frac_water and wcont_evap when there has been a change to wcont

	const double maxerr = 0.000001; // Max error allowed
															
	// Local storage
	double Faw_layer[NSOILLAYER];
	// available water for each soil layer (mm)
	double ice_layer[NSOILLAYER];
	// ice each soil layer (mm)

	Wtot = 0.0;

	for (int ly = IDX; ly < IDX + nsublayer1 + nsublayer2; ly++) {

		oob_check_wcont(wcont[ly - IDX]);

		if (wcont[ly - IDX] < 0.0 || wcont[ly - IDX] > 1.0) {
			dprintf("Soil::update_soil_water - bad wcont!\n");
			return;
		}

		// For the following, note that, for each layer: wcont = Faw_layer / soiltype.awc
		// Also, wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
		if (patch.stand.is_highlatitude_peatland_stand())
			Faw_layer[ly - IDX] = wcont[ly - IDX] * soiltype.awc_peat[ly - IDX]; // mm
		else
			Faw_layer[ly - IDX] = wcont[ly - IDX] * soiltype.awc[ly - IDX]; // mm

 		ice_layer[ly - IDX] = Frac_ice[ly] * Dz[ly]; // mm

		// Check balance
		double layerwater = Faw_layer[ly - IDX] + ice_layer[ly - IDX]; // mm

		if (aw_max[ly - IDX] - layerwater < -1.0 * maxerr) {
			dprintf("Soil::update_soil_water - error in a soil layer's water balance!\n");
			return;
		}

		// Now update Frac_water
		Frac_water[ly] = Faw_layer[ly - IDX] / Dz[ly];

		// For wetlands:
		if (ly < IDX + NACROTELM)
			Wtot += (Frac_water[ly] + Frac_ice[ly] + Fpwp_ref[ly]) * Dz[ly]; // Wtot includes ice

	} // for loop (ly)

	// Update wcont_evap
	// wcont_evap is the available, nonfrozen water in the top two layers (20 cm), divided by the capacity of those layers to hold water
	if (patch.stand.is_highlatitude_peatland_stand())
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc_peat[0] + soiltype.awc_peat[1]);
	else
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc[0] + soiltype.awc[1]);

	if (wcont_evap < 0.0 || wcont_evap > 1.0) {
		dprintf("Soil::update_soil_water - error in the wcont_evap calculation!\n");
		return;
	}

	oob_check_wcont(wcont_evap);

	return; // ...if no problems
}


/// Crank-Nicholson timestepper algorithm for temperature diffusion equation.
void cnstep_full(int layer0, double Di[NLAYERS], double dz[NLAYERS], double surf_temp, double dt,
	double pad_dz[PAD_LAYERS], double temp[NLAYERS], double pad_temp[PAD_LAYERS],
	double ki[NLAYERS], double ci[NLAYERS]) {

	// DESCRIPTION
	//     Crank-Nicholson timestepper for temperature diffusion equation.
	//     This routine performs a single timestep (of length dt) of the heat
	//     diffusion equation
	//
	//                        6 T    6  /      6 T \
	//                   c(z) --- = --- | k(z) --- |
	//                        6 t   6 z \      6 z /
	//
	//     where here, T(z, t) is the temperature, the '6's represent partial
	//     differentiation, and c(z)/k(z) is the (depth-dependent) heat capacities and thermal conductivities
	//     The boundary conditions are a prescribed surface
	//     temperature, and no heat flow at the bottom of the solution domain.
	//
	//     INPUT PARAMETERS:
	//
	//      NLAYERS      Total number of layers represented in temp, d and
	//                   dz arrays. Not all of the layers have to be active
	//                   at one time, but this parameter is used to constrain
	//                   the sizes of the depth-based arrays.
	//
	//      layer0       Index of the first active layer. The diffusion
	//                   equation will be solved for all layers in the range
	//                   layer0 - NLAYERS inclusive, and the surface temperature
	//                   boundary condition will be applied to the top of layer
	//                   index layer0.
	//
	//      temp         The temperatures in each layer at the start of the
	//                   timestep.
	//
	//      d            The thermal diffusivities in each layer.
	//
	//      k            The thermal conductivities in each layer.
	//
	//      c            The heat capacities in each layer.
	//
	//      dz           The thickness of each layer.
	//
	//      surf_temp    The surface temperature boundary condition.
	//
	//      dt           The timestep.
	//
	//      PAD_LAYERS   Number of layers to use below computation domain to
	//                   deal with bottom boundary condition.
	//
	//      pad_dz       Padding layer thicknesses (derived from calc_padding).
	//
	//      pad_temp     Temperature in padding layers.
	//
	//     NOTE ON UNITS: no units are specified here for any of these
	//     inputs; the only constraint on units is that the layer
	//     thicknesses, the diffusion constant and the timestep should be
	//     measured in consistent units. For instance, if the layer
	//     thicknesses are measured in metres, and the timestep in days, {
	//     the diffusion constants should have units of m2 / day.
	//
	//     OUTPUT PARAMETERS
	//
	//      temp         The temperatures in each layer at the end of the
	//                   timestep. Only layers layer0 - NLAYERS will have
	//                   valid temperature values after a to cnstep. All
	//                   layers with index less than layer0 will have missing
	//                   value temperatures set (a value of -1.0E35).
	//
	//      pad_temp     Temperature in padding layers at the end of the
	//                   timestep.

	// As cnstep, but now we feed in ki and ci

	// A fill-in for layers above layer0.
	double MISSING_VALUE;

	// Bedrock parameters from Chadburn et al. (2015)
	double Kbedrock = 8.6;			// W m-1 K-1
	double Cbedrock = 2100000.0;	// J m-3 K-1
	double Dbedrock = Kbedrock / Cbedrock * SECS_PER_DAY * MM2_PER_M2; // mm2 day-1

	// Layer counters: note that there are two different layer counting
	// schemes used, one for the input and output parameters (vectors of
	// length NLAYERS) and one for the values used in the Crank-Nicholson
	// solver (vectors of length active_layers)
	int layer, lidx;
	int active_layers;

	// Diffusion constants averaged over adjacent layers.
	double dplus;
	double dminus;

	// Layer-dependent weighting factors in Crank-Nicholson scheme.
	double dz_factor;
	double Cplus;
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double there, tminus, tplus;

	// Leading diagonal, left and right subdiagonals for Crank-Nicholson matrix.
	long double diag[active_layersmax];
	long double left[active_layersmax];
	long double right[active_layersmax];

	// Right hand side vector for Crank-Nicholson scheme equations.
	long double rhs[active_layersmax];

	// Solution vector for Crank-Nicholson scheme equations.
	long double solution[active_layersmax];

	// initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = 0.0;
	there = tminus = tplus = 0.0;

	for (int l = 0; l < active_layersmax; l++) {
		diag[l] = 0.0;
		left[l] = 0.0;
		right[l] = 0.0;
		rhs[l] = 0.0;
		solution[l] = 0.0;
	}

	// --- CODE STARTS HERE ---

	// Not all layers in the input arrays are actually used.  Calculate
	// how many layers we have to deal with here.
	MISSING_VALUE = -9999.0;

	active_layers = NLAYERS - layer0 + PAD_LAYERS;

	//   BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// End members for off-diagonal elements.
	left[0] = 0.0;
	right[active_layers - 1] = 0.0;

	// Process the active layers. 
	// The first time (lidx and layer = layer0) corresponds to the surface layer

	// lidx runs from 1 to active_layers == 
	for (lidx = 1; lidx <= active_layers; lidx++) {

		// Deal with different layer counting schemes.
		layer = lidx + layer0 - 1;
		// Minimum is layer0
		// Maximum is NLAYERS + PAD_LAYERS-1

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.  Diffusion
		// constants in the padding layers are calculated in a sensible way.

		// D+
		double small = 0.000001;

		// D+
		const int first_bedrock_padding_layer = PAD_LAYERS; // 1 bedrock layer

		if (layer == NLAYERS + PAD_LAYERS - 1)
			dplus = 0.0;				// BC2 - Bottom layer diffusion clamped to 0
		else if (layer >= NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer. 
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Dbedrock;			// The lowest padding layers are treated as bedrock
	 	else
			dplus = (1.0 / ci[layer]) * (ki[layer] * dz[layer] + ki[layer + 1] * dz[layer + 1]) / (dz[layer] + dz[layer + 1]);

		// D-

		if (layer == layer0)
			dminus = Di[layer];			// top layer
		else if (layer > NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer. 
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Dbedrock;			// The lowest padding layers are treated as bedrock
		else
			dminus = (1.0 / ci[layer]) * (ki[layer] * dz[layer] + ki[layer - 1] * dz[layer - 1]) / (dz[layer] + dz[layer - 1]);

		// Extract sensible values to use for temperature and diffusivity
		// of the current layer, the layer above and the layer below,
		// taking account of padding layers and end cases.

		// --- HERE ---
		if (layer < NLAYERS) {
			 // all soil layers
			dzhere = dz[layer];
			there = temp[layer];
		}
		else {
			// all padding layers
			dzhere = pad_dz[layer - NLAYERS];
			there = pad_temp[layer - NLAYERS];
		}

		// --- PLUS ---
		if (layer < NLAYERS - 1) {
			// all soil layers apart from the bottom soil layer
			dzplus = dz[layer + 1];
			tplus = temp[layer + 1];
		}
		else if (layer < NLAYERS + PAD_LAYERS - 1) {
			// bottom soil layer and all 
			// padding layers apart from the very bottom padding layer 
			dzplus = pad_dz[layer - NLAYERS + 1];
			tplus = pad_temp[layer - NLAYERS + 1];
		}
		else {
			// bottom padding layer
			dzplus = pad_dz[PAD_LAYERS - 1];
			tplus = pad_temp[PAD_LAYERS - 1];
		}

		// --- MINUS ---
		if (layer == layer0) {
			// top layer
			dzminus = dz[layer];
			tminus = temp[layer];
		}
		else if (layer <= NLAYERS) {
			// all soil layers and the first padding layer
			dzminus = dz[layer - 1];
			tminus = temp[layer - 1];
		}
		else {
			// all padding layers apart from the top padding layer
			dzminus = pad_dz[layer - NLAYERS - 1];
			tminus = pad_temp[layer - NLAYERS - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.

		// DIAG
		if (lidx == 1)
			diag[0] = 1.0; // BC1 - top layer should be (1,0,...,0)
		else
			diag[lidx - 1] = 1.0 + Cplus + Cp_minus;

		// LEFT & RIGHT
		if (lidx < active_layers) {

			if (lidx > 1)
				right[lidx - 1] = -Cplus;
			else
				right[lidx - 1] = 0.0; // i.e. BC1, where top layer == (1,0,..,0)

			// left[0] is set above.
			if (lidx > 1)
				left[lidx - 1] = -Cp_minus;
		}

		if (lidx == active_layers) left[lidx - 1] = -Cp_minus;

		// RHS
		// Calculate right hand side vector values.
		if (lidx == 1)
			rhs[0] = surf_temp;
		else if (lidx == active_layers) // Cplus == 0 here anyway
			rhs[lidx - 1] = (1.0 - Cp_minus) * there + Cp_minus * tminus;
		else
			rhs[lidx - 1] = (1.0 - Cplus - Cp_minus) * there +
			Cplus * tplus + Cp_minus * tminus;
	} // end for

	//   SOLVE TRIDIAGONAL SYSTEM

	tridiag(active_layers, left, diag, right, rhs, solution);

	// Numerical tests

	if (DEBUG_SOIL_TEMPERATURE) {

		int testrow = 10;
		double checksum = left[testrow] * solution[testrow - 1] + diag[testrow] * solution[testrow] +
			right[testrow] * solution[testrow + 1] - rhs[testrow];

		if (fabs(checksum) > 0.0000000001)
			dprintf("%s\n", "Bad checksum after tridiag - test1");
	}

	//   FORMAT OUTPUT PARAMETERS

	// Transfer the solution to the temperature array.  Also update
	// the padding layer temperatures so they can be used for the
	// next timestep.
	for (lidx = 0; lidx<NLAYERS; lidx++) {

		// remove temps < 1e-12.
		if (verysmall(solution[lidx])) solution[lidx] = 0.0;

		if (lidx < layer0)
			temp[lidx] = MISSING_VALUE;
		else
			temp[lidx] = solution[lidx - layer0];
	}

	for (lidx = 0; lidx<PAD_LAYERS; lidx++) {
		pad_temp[lidx] = solution[lidx + NLAYERS - layer0];
	}
}

/// Crank-Nicholson timestepper algorithm for temperature diffusion equation.
void cnstep(int layer0, double Di[NLAYERS], double dz[NLAYERS], double surf_temp, double dt,
	double pad_dz[PAD_LAYERS], double temp[NLAYERS], double pad_temp[PAD_LAYERS]) {
	// DESCRIPTION
	// Crank-Nicholson timestepper for temperature diffusion equation.
	// This routine performs a single timestep (of length dt) of the heat
	// diffusion equation
	//
	//                        6 T    6  /      6 T \
	//                        --- = --- | D(z) --- |
	// 						  6 t   6 z \      6 z /
	//
	// where here, T(z, t) is the temperature, the '6's represent partial
	// differentiation, and D(z) is the (depth-dependent) diffusion
	// constant.  The boundary conditions are a prescribed surface
	// temperature, and no heat flow at the bottom of the solution domain.
	//
	// INPUT PARAMETERS:
	//
	// NLAYERS      Total number of layers represented in temp, d and
	//              dz arrays.  Not all of the layers have to be active
	//              at one time, but this parameter is used to constrain
	//              the sizes of the depth-based arrays.
	//
	// layer0       Index of the first active layer.  The diffusion
	//              equation will be solved for all layers in the range
	//              layer0 - NLAYERS inclusive, and the surface temperature
	//              boundary condition will be applied to the top of layer
	//              index layer0.
	//
	// temp         The temperatures in each layer at the start of the
	//              timestep.
	//
	// Di           The thermal diffusivities in each layer. (m2/day)
	//
	// dz           The thickness of each layer.
	//
	// surf_temp    The surface temperature boundary condition.
	//
	// dt           The timestep.
	//
	// PAD_LAYERS   Number of layers to use below computation domain to
	//              deal with bottom boundary condition.
	//
	// pad_dz       Padding layer thicknesses (derived from calc_padding).
	//
	// pad_temp     Temperature in padding layers.
	//
	// NOTE ON UNITS: no units are specified here for any of these
	// inputs; the only constraint on units is that the layer
	// thicknesses, the diffusion constant and the timestep should be
	// measured in consistent units.  For instance, if the layer
	// thicknesses are measured in metres, and the timestep in days, {
	// the diffusion constants should have units of m2 / day.
	//
	// OUTPUT PARAMETERS
	//
	// temp         The temperatures in each layer at the end of the
	//              timestep.  Only layers layer0 - NLAYERS will have
	//              valid temperature values after a  to cnstep.  All
	//              layers with index less than layer0 will have missing
	//              value temperatures set (a value of -1.0E35).
	//
	// pad_temp     Temperature in padding layers at the end of the
	//              timestep.


	// A fill-in for layers above layer0.
	double MISSING_VALUE;

	// Bedrock parameters from Chadburn et al. (2015)
	double Kbedrock = 8.6;			// W m-1 K-1
	double Cbedrock = 2100000.0;	// J m-3 K-1
	double Dbedrock = Kbedrock / Cbedrock * SECS_PER_DAY * MM2_PER_M2; // mm2 day-1

	// Layer counters: note that there are two different layer counting
	// schemes used, one for the input and output parameters (vectors of
	// length NLAYERS) and one for the values used in the Crank-Nicholson
	// solver (vectors of length active_layers)
	int layer, lidx;
	int active_layers;

	// Diffusion constants averaged over adjacent layers.
	double dplus;
	double dminus;

	// Layer-dependent weighting factors in Crank-Nicholson scheme.
	double dz_factor;
	double Cplus;
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double there, tminus, tplus;

	// Leading diagonal, left and right subdiagonals for Crank-Nicholson
	// matrix.
	long double diag[active_layersmax];
	long double left[active_layersmax];
	long double right[active_layersmax];

	// Right hand side vector for Crank-Nicholson scheme equations.
	long double rhs[active_layersmax];

	// Solution vector for Crank-Nicholson scheme equations.
	long double solution[active_layersmax];

	// initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = 0.0;
	there = tminus = tplus = 0.0;

	// --- CODE STARTS HERE ---

	// Not all layers in the input arrays are actually used.  Calculate
	// how many layers we have to deal with here.
	MISSING_VALUE = -9999.0;

	active_layers = NLAYERS - layer0 + PAD_LAYERS;

	//   BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// End members for off-diagonal elements.
	left[0] = 0.0;
	right[active_layers - 1] = 0.0;

	// Process the active layers. 
	// The first time (lidx and layer = layer0) corresponds to the surface layer

	for (lidx = 1; lidx <= active_layers; lidx++) {

		if (lidx == 1) { // moved from above
			diag[lidx-1] = 0.0;
			left[lidx-1] = 0.0;
			right[lidx-1] = 0.0;
			rhs[lidx-1] = 0.0;
			solution[lidx-1] = 0.0;
		}

		// Deal with different layer counting schemes.
		layer = lidx + layer0 - 1;
		// Minimum is layer0
		// Maximum is NLAYERS + PAD_LAYERS-1

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.  Diffusion
		// constants in the padding layers are calculated in a sensible way.

		// D+
		const int first_bedrock_padding_layer = 3; // 3 layers

		if (layer == NLAYERS + PAD_LAYERS - 1)
			dplus = 0.0;				// BC2 - Bottom layer diffusion clamped to 0
		else if (layer >= NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer. 
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Dbedrock;			// The lowest padding layers are treated as bedrock
		else
			dplus = 0.5 * (Di[layer] + Di[layer + 1]);

		// D-

		if (layer == layer0)
			dminus = Di[layer]; // top layer
		else if (layer > NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer.
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Dbedrock;			// The lowest padding layers are treated as bedrock
		else
			dminus = 0.5 * (Di[layer] + Di[layer - 1]); // soil layers

		// Extract sensible values to use for temperature and diffusivity
		// of the current layer, the layer above and the layer below,
		// taking account of padding layers and end cases.

		// --- HERE ---
		if (layer < NLAYERS) {
			// all soil layers
			dzhere = dz[layer];
			there = temp[layer];
		}
		else {
			// all padding layers
			dzhere = pad_dz[layer - NLAYERS];
			there = pad_temp[layer - NLAYERS];
		}

		// --- PLUS ---
		if (layer < NLAYERS - 1) {
			// all soil layers apart from the bottom soil layer
			dzplus = dz[layer + 1];
			tplus = temp[layer + 1];
		}
		else if (layer < NLAYERS + PAD_LAYERS - 1) {
			// bottom soil layer and all 
			// padding layers apart from the very bottom padding layer 
			dzplus = pad_dz[layer - NLAYERS + 1];
			tplus = pad_temp[layer - NLAYERS + 1];
		}
		else {
			// bottom padding layer
			dzplus = pad_dz[PAD_LAYERS - 1];
			tplus = pad_temp[PAD_LAYERS - 1];
		}

		// --- MINUS ---
		if (layer == layer0) {
			// top layer
			dzminus = dz[layer];
			tminus = temp[layer];
		}
		else if (layer <= NLAYERS) {
			// all soil layers and the first padding layer
			dzminus = dz[layer - 1];
			tminus = temp[layer - 1];
		}
		else {
			// all padding layers apart from the top padding layer
			dzminus = pad_dz[layer - NLAYERS - 1];
			tminus = pad_temp[layer - NLAYERS - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.

		// DIAG
		if (lidx == 1)
			diag[0] = 1.0; // BC1 - top layer should be (1,0,...,0)
		else
			diag[lidx - 1] = 1.0 + Cplus + Cp_minus;

		// LEFT & RIGHT
		if (lidx < active_layers) {

			if (lidx > 1)
				right[lidx - 1] = -Cplus;
			else
				right[lidx - 1] = 0.0; // i.e. BC1, where top layer == (1,0,..,0)

			// left[0] is set above.
			if (lidx > 1)
				left[lidx - 1] = -Cp_minus;
		}

		if (lidx == active_layers) left[lidx - 1] = -Cp_minus;

		// RHS
		// Calculate right hand side vector values.
		if (lidx == 1)
			rhs[0] = surf_temp;
		else if (lidx == active_layers) // Cplus == 0 here anyway
			rhs[lidx - 1] = (1.0 - Cp_minus) * there + Cp_minus * tminus;
		else
			rhs[lidx - 1] = (1.0 - Cplus - Cp_minus) * there +
			Cplus * tplus + Cp_minus * tminus;
	} // end for

	// SOLVE TRIDIAGONAL SYSTEM

	tridiag(active_layers, left, diag, right, rhs, solution);

	// Numerical test

	if (DEBUG_SOIL_TEMPERATURE) {

		int testrow = 10;
		double checksum = left[testrow] * solution[testrow - 1] + diag[testrow] * solution[testrow] +
			right[testrow] * solution[testrow + 1] - rhs[testrow];

		if (fabs(checksum) > 0.0000000001)
			dprintf("%s\n", "Bad checksum after tridiag - test1");
	}

	//   FORMAT OUTPUT PARAMETERS

	// Transfer the solution to the temperature array.  Also update
	// the padding layer temperatures so they can be used for the
	// next timestep.
	for (lidx = 0; lidx<NLAYERS; lidx++) {

		if (lidx < layer0)
			temp[lidx] = MISSING_VALUE;
		else
			temp[lidx] = solution[lidx - layer0];

		if (lidx < PAD_LAYERS) // was in the loop below
			pad_temp[lidx] = solution[lidx + NLAYERS - layer0];
	}
}


void calc_padding(double dz, double pad_dz[PAD_LAYERS]) {

	// DESCRIPTION
	// Calculate padding layer depths for Crank-Nicholson bottom boundary
	// condition setup. The padding layer thicknesses are calculated as
	// a geometric sequence with a given total padding thickness and
	// number of layers. This gives a smooth increase in layer depth
	// away from the calculation domain.

	// INPUT PARAMETERS
	// dz;					// Thickness of last double layer.
	// pad_dz[PAD_LAYERS]	// Output padding layer thicknesses [mm].

	// Function, derivative and parameter for geometric factor
	// calculation.
	// double func, dfunc, k

	// use the global variables PAD_DEPTH and PAD_LAYERS
	double pad_depth = PAD_DEPTH;

	// Iterative improvements of geometric factor
	double k_guess, k_last;

	// Layer counter.
	int layer;

	// Statement functions for function and derivative for Newton
	// solver for geometric sequence factor.

	// Initial state - don't choose 1.0 as the initial value for
	// k_guess, since it's a trivial solution//
	k_guess = 10.0;
	k_last = 0.0;

	double func, dfunc;

	// Do Newton-Raphson solve to find geometric sequence factor.
	while (fabs(k_guess - k_last) > 1E-2) {
		k_last = k_guess;

		func = pow(k_guess, (PAD_LAYERS + 1)) - k_guess * (1.0 + pad_depth / dz) + pad_depth / dz;
		dfunc = (PAD_LAYERS + 1) * pow(k_guess, PAD_LAYERS) - (1.0 + pad_depth / dz);
		k_guess = k_guess - func / dfunc;
	}

	// Build geometric sequence of padding layer depths.
	for (layer = 0; layer<PAD_LAYERS; layer++) {
		pad_dz[layer] = dz * pow(k_guess, layer + 1);
	}
}


void Soil::init_hydrology_variables() {

	// DESCRIPTION
	// Initialisation of hydrology. 
	// Called ONCE, at the beginning of every stand's simulation
	// The initialisation relies on the following definitions:

	// sub_water[NSUBLAYERS_ACRO] - the volumetric water content in the NSUBLAYERS_ACRO of the acrotelm
	// whc[NSOILLAYER] - available water holding capacity of soil layers [0=upper layer] [mm], taking into
	// account the unavailability of frozen water. Default value: soiltype.awc[]
	// aw_max[NSOILLAYER] - Max water (mm) that can be held in each layer
	// alwhc[NLAYERS] - the volumetric liquid water content. A fraction. 
	// Considers the entire (awc + Fpwp) volumetric water content MINUS the ice fraction. Updated daily.
	// alwhc_init[NLAYERS] - the initial volumetric liquid water content. A fraction. Considers the entire
	// (awc + Fpwp) volumetric water content MINUS the ice fraction.

	// Procedure:
	// Initialise whc[] and the total volumetric liquid water content, alwhc[]
	// Save the initial values - alwhc_init 
	// Initialise aw_max[], the max water (mm) that can be held in each layer

	for (int ly = 0; ly < nsublayer1 + nsublayer2; ly++) {
		
		if (patch.stand.is_highlatitude_peatland_stand()) {

			// Override hydrology properties for peatland cells
			if (ly < NACROTELM) {
				alwhc[ly] = acro_por - peat_wp;
				soiltype.awc_peat[ly] = Dz_acro * (acro_por-peat_wp); // mm
				soiltype.wsats_peat[ly] = acro_por * Dz_acro;
			}
			else {
				alwhc[ly] = cato_por - peat_wp;
				soiltype.awc_peat[ly] = Dz_cato * (cato_por-peat_wp); // mm
				soiltype.wsats_peat[ly] = cato_por * Dz_cato;
			}
			
			whc[ly] = soiltype.awc_peat[ly]; // mm
		}
		else {

			// Not a high latitude peatland
			if (ly < nsublayer1)
				alwhc[ly] = soiltype.awc[ly] / (SOILDEPTH_UPPER / NSOILLAYER_UPPER);
			else
				alwhc[ly] = soiltype.awc[ly] / (SOILDEPTH_LOWER / NSOILLAYER_LOWER);

			whc[ly] = soiltype.awc[ly]; // mm
		}
		
		alwhc_init[ly] = alwhc[ly];
		
		// Max water (mm) that can be held in each layer, above the pwp
		aw_max[ly] = alwhc_init[ly] * Dz[ly+IDX]; // mm
	}
}


void Soil::update_acrotelm_co2(double atmo_co2) {

	// DESCRIPTION
	// Update the CO2 concentration in the acrotelm water, after Wania et al. (2009)
	// Varies between atmospheric CO2 level when awtp = -300mm, to PORE_WATER_CO2 when awtp = 0, 
	// and greater if there's standing water.
	// Called once daily for wetlands, after the water table has been updated

	// LPJ-WHy comments from Rita Wania:
	/*
	     LPJ-WHy uses a different atmospheric CO2 concentration for mosses
	     as they can access the CO2 dissolved in acrotelm water. If the
	     water table is high, mosses can access all of the acrotelm CO2,
	     but as the water table drops the CO2 concentrations available to
	     mosses is a mixture between atmospheric CO2 and acrotelm CO2.
	     NOTE: for ideal gases, ppmv = not just mimol L-1 but also
	     mimol mol-1.
	     The mean annual wtd of the previous year is used here, as the
	     water table is only calculated after the co2 concentration is
	     needed
	*/

	acro_co2 = min(PORE_WATER_CO2, PORE_WATER_CO2 + ((PORE_WATER_CO2 - atmo_co2) / acro_depth) * awtp);
}


bool Soil::update_layer_water_content(int day) {

	// DESCRIPTION
	// Private helper method to update the soil water contents in the soil layers. 
	// Called each day.

	// whc[] is water above the pwp which hasn't frozen [mm]
	// Its maximum value is awc[] (mm), i.e. whc[later] <= awc[layer], for all layers

	const double maxerr = 0.000001; // Max error allowed

	// Local storage
	double Faw_layer[NSOILLAYER];
	// available water for each soil layer (mm)
	double ice_layer[NSOILLAYER];
	// mm of ice in each soil layer (mm)

	// Peatland total soil water content, incl ice and water below the pwp
	Wtot = 0.0; // mm

	for (int ly = IDX; ly < IDX + nsublayer1 + nsublayer2; ly++) {

		// Update the total available volumetric water content (mm) of the layers
		Faw_layer[ly-IDX] = Frac_water[ly] * Dz[ly]; // mm
		ice_layer[ly-IDX] = Frac_ice[ly] * Dz[ly]; // mm

		// Check balance
		double layerwater = Faw_layer[ly-IDX] + ice_layer[ly-IDX]; // mm

		if (aw_max[ly-IDX] - layerwater < -1.0 * maxerr) {
			dprintf("Soil::update_layer_water_content - error in a soil layer's water balance!\n");
			return false;
		}

		// Subtract the ice/frozen fraction from the initial value of the capacity for this layer to hold water
		alwhc[ly-IDX] = alwhc_init[ly-IDX] - Frac_ice[ly]; 

		// For wetlands:
		// total ice in layer = Frac_ice[i] + Fpwp_ref[i] - Frac_water_belowpwp[i];
		// total water in layer = Frac_water[i] + Frac_water_belowpwp[i];
		// So, the total ice + water = Frac_ice[i] + Fpwp_ref[i] + Frac_water[i] 
		if (ly < IDX + NACROTELM)
			Wtot += (Frac_water[ly] + Frac_ice[ly] + Fpwp_ref[ly]) * Dz[ly]; // Wtot includes ice

		if (patch.stand.is_highlatitude_peatland_stand())
			whc[ly-IDX] = max(0.0, aw_max[ly-IDX] - ice_layer[ly-IDX]); // mm
		else
			whc[ly-IDX] = max(0.0, soiltype.awc[ly-IDX] - ice_layer[ly-IDX]); // mm 
		
		// Avoid dividing by small numbers later
		if (fabs(whc[ly-IDX]) < maxerr)
			whc[ly-IDX] = 0.0;

		// wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
		if (patch.stand.is_highlatitude_peatland_stand())
			wcont[ly-IDX] = Faw_layer[ly-IDX] / soiltype.awc_peat[ly-IDX];
		else
			wcont[ly-IDX] = Faw_layer[ly-IDX] / soiltype.awc[ly-IDX];

		oob_check_wcont(wcont[ly - IDX]);

		if (wcont[ly-IDX] < 0.0 || wcont[ly-IDX] > 1.0) {
			dprintf("Soil::update_layer_water_content - error in the wcont calculation!\n");
			return false;
		}
	} // for loop (ly)

	// wcont_evap is the available, nonfrozen water in the tope two layers (20 cm), divided by the capacity of those layers to hold water
	if (patch.stand.is_highlatitude_peatland_stand())
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc_peat[0] + soiltype.awc_peat[1]);
	else 
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc[0] + soiltype.awc[1]);
	
	oob_check_wcont(wcont_evap);

	if (wcont_evap < 0.0 || wcont_evap > 1.0) {
		dprintf("Soil::update_layer_water_content - error in the wcont_evap calculation!\n");
		return false;
	}

	return true; // ...if no problems
}


// Private helper function
void Soil::update_from_yesterday() {

	// DESCRIPTION
	// Update T_soil and Frac_ice using this function. 
	// More variables could be added here.

	// Frac_ice_yesterday and T_soil_yesterday record yesterday's values, not just those on Dec 31.
	for (int i = 0; i<NLAYERS; i++) {
		T_soil[i] = T_soil_yesterday[i];
		Frac_ice[i] = Frac_ice_yesterday[i];
	}
}


void Soil::update_snow_properties(const int& daynum, const double& dailyairtemp, double& Dsnow, double& Csnow, double& Ksnow) {

	// DESCRIPTION
	// Hillel (1982) gives values for porosity, Csnow and Ksnow.       
	// Here we update snow density and use that to estimate Csnow and Ksnow 

	// Csnow:		heat capacity of snow [J m-3 K-1]
	// Ksnow:		thermal conductivity of snow [W m-1 K-1]
	// Dsnow:		thermal conductivity of snow [mm2 day-1]

	bool snow_active_old = snow_active;

	// Are there snow layers?
	if (snowpack > 1) // mm water
		snow_active = true;
	else
		snow_active = false;

	// No snow, so return
	if (!snow_active) {
		snow_active_layers = 0;
		return;
	}

	// new snow?
	if (!snow_active_old)
		snowdens = snowdens_start;

	// Update snowdens using the Wania et al. simple snow compaction scheme

	// We assume fresh snow unless snow_days > 0.75 * 365;
	if (date.year >= 1 && snow_days >= 0.75 * snow_days_prev) {
		snowdens = snowdens_start + (snowdens_end - snowdens_start) / (snow_days_prev * 0.25) * (snow_days - 0.75 * snow_days_prev);
		 
		// Restrict density to the range [snowdens_end,snowdens_start]
		snowdens = max(min(snowdens, snowdens_end), snowdens_start);
	}
	else
		snowdens = snowdens_start;

	// Override the above with a fixed snow density?
	if (snowdensityconstant) {
		snowdens = 250.0; // kg m-3, JSBACH value from Ekici et al (2015)
	}

	// Note that K/C = D has units mm2/day, as above
	// Wania values
	// Csnow from Fukusako, Eq. 2. It's the same as for ice
	static const double ice_density = 917.0; // kg m-3
	Csnow = (0.185 + 0.689 * (K2degC + dailyairtemp) * 0.01) * J_PER_KJ; // J kg-1 K-1
	Csnow *= ice_density; // conversion to volumetric heat capacity (J m-3 K-1) - should be 1,900,000 J m-3 K-1 approx.

	double snowdens_scaled = snowdens / rho_H2O;
	// From Sturm et al. 1997, and in Ling and Zhang, 2006
	if (snowdens > 156.0)
		Ksnow = 0.138 - 1.01 * snowdens_scaled + 3.233 * snowdens_scaled * snowdens_scaled;
	else
		Ksnow = 0.023 + 0.234 * snowdens_scaled;
	 
	Dsnow = Ksnow / Csnow * SECS_PER_DAY * MM2_PER_M2; // mm2 day-1

	// Unit conversions - needed in cnstep_full:
	// Conversion from (W m-1 K-1) TO (J day-1 mm-1 K-1) 
	Ksnow *= SECS_PER_DAY * M_PER_MM;
	// Conversion from (J m-3 K-1) TO (J mm-3 K-1)
	Csnow *= M3_PER_MM3;
}


// A valid layer number?
bool Soil::valid_layer_num(const int& ngr) const {

	if (ngr != (SOILDEPTH_UPPER + SOILDEPTH_LOWER) / Dz_soil)
		return false;
	else
		return true;
}


void Soil::update_soil_diffusivities(const int& daynum, bool ansoln) {

	// DESCRIPTION
	// Called daily to update diffusivities, heat capacities and thermal conductivities of the 
	// soil layers below the mixed/moss and snow layers 

	double lKfwater;					// log of therm. cond. for the water fraction
	double lKfice;						// log of therm. cond. for the ice fraction
	double lKforg;						// log of therm. cond. for the org. fraction
	double lKfpeat;						// log of therm. cond. for the peat fraction
	double lKfmin;						// log of therm. cond. for the min. fraction
	double ThermCond_total[NLAYERS];	// total thermal condunctivity
	double HeatCapacity_total[NLAYERS];	// total heat capacity
	double Ftot;						// total of solid fraction in soil

	for (int i = IDX; i<NLAYERS; i++) {

		if (ansoln) Frac_water[i] = 0.0; // No water if analytic soln. sought

		// take water and ice below pwp into account
		double totaliceinlayer = Frac_ice[i] + Fpwp_ref[i] - Frac_water_belowpwp[i];
		double totalwaterinlayer = Frac_water[i] + Frac_water_belowpwp[i];

		Frac_air[i] = por[i] - totalwaterinlayer - totaliceinlayer;

		if (Frac_air[i] < -0.00001)
			fail("Illegal Frac_air value in Soil::update_soil_diffusivities: %g", Frac_air[i]);

		// Catch rounding errors for Frac_air
		if (Frac_air[i] < 0.0 && Frac_air[i] > -0.00001)
			Frac_air[i] = 0.0;

		if (ansoln) Frac_air[i] = 0.0; // No air if analytic soln. sought

		HeatCapacity_total[i] = Frac_min[i] * Cp_min + Frac_org[i] * Cp_org + totaliceinlayer * Cp_ice +
			 totalwaterinlayer * Cp_water + Frac_peat[i] * Cp_peat + Frac_air[i] * Cp_air;

		// Update the THERMAL CONDUCTIVITY K - which follows
		// Granberg et al. 1999, who refers to Farouki, 1986

		Ftot = totalwaterinlayer + totaliceinlayer + Frac_org[i] + Frac_min[i] + Frac_peat[i]; 

		if (Ftot + Frac_air[i] > 1.000001)
			fail("Illegal Ftot (%g) + Frac_air (%g) value in Soil::update_soil_diffusivities: %g", Ftot, Frac_air[i], Ftot +Frac_air[i]);

		lKforg = (Frac_org[i] / Ftot) * lKorg;
		lKfpeat = (Frac_peat[i] / Ftot) * lKpeat;
		lKfmin = (Frac_min[i] / Ftot) * lKmin;
		lKfwater = (totalwaterinlayer / Ftot) * lKwater;
		lKfice = (totaliceinlayer / Ftot) * lKice;

		ThermCond_total[i] = Frac_air[i] * Kair + Ftot * exp(lKfwater + lKfice + lKforg + lKfmin + lKfpeat);

		// Update THERMAL DIFFUSIVITY for each SOIL layer
		// multiply by 86400 to get from m2 s-1 to m2 d-1
		// multiply by 1E6 to get from m2 d-1 to mm2 d-1
		Di[i] = ThermCond_total[i] / HeatCapacity_total[i] * SECS_PER_DAY * MM2_PER_M2;	

		// Conversion from (W m-1 K-1) TO (J day-1 mm-1 K-1)
		Ki[i] = ThermCond_total[i] * SECS_PER_DAY * M_PER_MM;
		// Conversion from (J m-3  K-1) TO (J mm-3 K-1)
		Ci[i] = HeatCapacity_total[i] * M3_PER_MM3;

		// Note that K/C = D has units mm2/day, as above
	} // i 
}


void Soil::update_ice_fraction(const int& daynum, const int& MIDX) {

	// DESCRIPTION
	// Calculate WATER AND ICE fraction in each layer

	// To know how much ice is formed or melted at each time step,
	// The temperature change from the previous to this time step is first 
	// calculated. If not all of the available water is frozen, the
	// temperature will be held at 0 deg until everything is frozen (and
	// vice versa for the melting). The amount of heat released when 
	// freezing or absorbed when thawing is to be found. To do this
	// the total heat capacity of the soil is multiplied with delta_T. Since
	// we are working in m-3 here, we don't need to multiply by the volume
	// or mass, because we want to get the fraction of water freezing or
	// ice thawing, not the absolute value. 

	// Determine where to start:
	// if there's standing water, include that layer for the phase change 
	// calculation, if not, start with the top soil layer


	double Ffreez;	// fraction of water which freezes
	double Fthaw;	// fraction of ice which thaws
	double energy;	// Energy per volume [J m-3]

	const double SMALLVOLFRAC = 0.0000001;	// Minimum volume fraction [m-3] allowed

	double delta_T;	// temperature difference (current - previous)

	for (int i = MIDX; i<NLAYERS; i++) {

		if (firstTempCalc)
			delta_T = 0.0;

		double heat_capacity_layer = Ci[i] * MM3_PER_M3; // J m-3 K-1

		double waterabovepwp = Frac_ice[i] + Frac_water[i]; 
		double ice_thislayer = Frac_ice[i] + (Fpwp_ref[i]-Frac_water_belowpwp[i]);
		double water_thislayer = Frac_water[i] + Frac_water_belowpwp[i];
		// The total amount of ice and water should not change in this routine
		double wtotal_initial = ice_thislayer + water_thislayer;

		// Add the water (not ice) below wp before calculating phase change
		Frac_water[i] += Frac_water_belowpwp[i];
		Frac_ice[i] += (Fpwp_ref[i]-Frac_water_belowpwp[i]);

		// Initialise temperatures on first day
		if (daynum == 0 && date.year == FIRST_FREEZE_YEAR)
			T_old[i] = T_soil[i];

		delta_T = T_soil[i] - T_old[i];

		// FREEZING
		if (T_soil[i] < 0.0 && T_soil[i] < T_old[i]) { // Wania conditions 
			Fthaw = 0.0;
			if (Frac_water[i] > 0.0 /*&& T_old[i] < 1.0*/) {
				energy = Cp_water * fabs(delta_T); //energy in J m-3 
				Ffreez = energy / Lheat; //unitless fraction

				if (Frac_water[i] > Ffreez) {
					Frac_ice[i] += Ffreez;
					Frac_water[i] -= Ffreez;

					// How much of the water is below pwp?
					if (Frac_water[i] < Fpwp_ref[i]) {
						Frac_water_belowpwp[i] = Frac_water[i];  
						Frac_water[i] = 0.0; // So there is liquid water, but it's all below the pwp
					}
					else {
						Frac_water_belowpwp[i] = Fpwp_ref[i];  
						Frac_water[i] -= Fpwp_ref[i]; // So there is liquid water above the pwp
					}

					T_soil[i] = 0.0;

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_water[i] < SMALLVOLFRAC)
							Frac_water[i] = 0.0;
						if (Frac_ice[i] < SMALLVOLFRAC) 
							Frac_ice[i] = 0.0;
						if (Frac_water_belowpwp[i] < SMALLVOLFRAC) 
							Frac_water_belowpwp[i] = 0.0;
					}
				}
				else {

					// It is cold enough to freeze all the water and decrease the temperature below 0, 
					// since the rest of the water gets frozen now. However, it should never result in 
					// temperatures colder than the values coming from the Crank-Nicholson scheme.

					Frac_ice[i] += Frac_water[i];
					
					// Update soil temperature in this layer
					T_soil[i] = -(Ffreez - Frac_water[i]) * Lheat / Cp_water;

					// Note that this gives the same result as the following method:
					// double temp_rise = Frac_water[i] * Lheat / Cp_water; // temp rise resulting from freezing water.
					// T_soil[i] += temp_rise;
					
					T_soil[i] = min(0.0, T_soil[i]); // Ensures we don't go back over 0 degrees

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_ice[i] < SMALLVOLFRAC)
							Frac_ice[i] = 0.0;
						if (Frac_water_belowpwp[i] < SMALLVOLFRAC) 
							Frac_water_belowpwp[i] = 0.0;
					}

					Frac_water[i] = 0.0;
					Frac_water_belowpwp[i] = 0.0; // All water is frozen, even that below the pwp
				}
			}
			else {
				Ffreez = 0.0;
				Frac_water[i] -= Frac_water_belowpwp[i];	
			}	
		} // end of FREEZING

		// THAWING
		else if (T_soil[i] >= 0.0  && T_soil[i] > T_old[i]) { // Wania conditions
			Ffreez = 0.0;
			if (Frac_ice[i] > 0.0) {
				energy = Cp_water * fabs(delta_T); //energy in J m-3
				Fthaw = energy / Lheat; //unitless fraction
				if (Frac_ice[i] > Fthaw) {

					// There is some ice left, even after melt
					Frac_water[i] += Fthaw;

					// How much of the water is below the pwp now?
					if (Frac_water[i] < Fpwp_ref[i]) {
						Frac_water_belowpwp[i] = Frac_water[i];  
						Frac_water[i] = 0.0; // So there is liquid water, but it's all below the pwp
					}
					else {
						Frac_water_belowpwp[i] = Fpwp_ref[i];  
						Frac_water[i] -= Fpwp_ref[i]; // So there is liquid water above the pwp
					}

					Frac_ice[i] -= Fthaw;
					T_soil[i] = 0.0;

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_ice[i] < SMALLVOLFRAC) 
							Frac_ice[i] = 0.0;
						if (Frac_water[i] < SMALLVOLFRAC) 
							Frac_water[i] = 0.0;
						if (Frac_water_belowpwp[i] < SMALLVOLFRAC) 
							Frac_water_belowpwp[i] = 0.0;
					}
				}
				else {
					// There is more energy available than required for
					// melting the rest of the ice. Therefore the
					// left-over energy will increase the temperature.

					Frac_water[i] += Frac_ice[i];
					Frac_water_belowpwp[i] = Fpwp_ref[i]; // All water is unfrozen
					Frac_water[i] -= Frac_water_belowpwp[i]; // Only store liquid water above the pwp

					// Update soil temperature in this layer
					T_soil[i] = (Fthaw - Frac_ice[i]) * Lheat / Cp_water;					

					// Note that this gives the same result as the following method:
					// double temp_reduction = Frac_ice[i] * Lheat / Cp_water; // temp reduction resulting from thawing ice.
					// T_soil[i] -= temp_reduction;

					T_soil[i] = max(0.0, T_soil[i]); // Ensures we don't go back under 0 degrees

					Frac_ice[i] = 0.0;

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_water[i] < SMALLVOLFRAC) 
							Frac_water[i] = 0.0;
					}
				}
			}
			else {
				Fthaw = 0.0;
				Frac_ice[i] = 0.0;
				Frac_water[i] -= Frac_water_belowpwp[i];	
			}
		} // end of THAWING
		else {
			// Neither freezing nor thawing, so revert to Frac_water
			Frac_water[i] -= Frac_water_belowpwp[i];	
		}

		// Subtract ice below pwp from Frac_ice 
		Frac_ice[i] -= (Fpwp_ref[i]-Frac_water_belowpwp[i]);

		if (Frac_ice[i] < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid ice volume content %g\n", Frac_ice[i]);

		if (Frac_water[i] < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid volumetric water content %g\n", Frac_water[i]);

		// Remove tiny amounts of ice and water for numerical stability
		if (DEBUG_SOIL_WATER) {
			// Tidy up due to rounding errors
			if (Frac_ice[i] < 0.0 && Frac_ice[i] > -1.0 * SMALLVOLFRAC)
				Frac_ice[i] = 0.0;

			if (Frac_water[i] < 0.0 && Frac_water[i] > -1.0 * SMALLVOLFRAC)
				Frac_water[i] = 0.0;
		}

		// water balance checks:
		double ice_thislayer_new = Frac_ice[i] + (Fpwp_ref[i]-Frac_water_belowpwp[i]);
		double water_thislayer_new = Frac_water[i] + Frac_water_belowpwp[i];
		double wtotal_final = ice_thislayer_new + water_thislayer_new;

		if ((fabs(wtotal_initial - wtotal_final)) > SMALLVOLFRAC)
			fail("Water mass balance during freezing or thawing in Soil::update_ice_fraction. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);

		// Update T_old after possible phase change
		T_old[i] = T_soil[i];

	} // end of loop through layers
}


void Soil::update_layer_fractions(const int& daynum, const int& mixedl, const int& MIDX) {

	// DESCRIPTION
	// Inializes and updates fractions of water, mineral content, peat content, organic content, as well as porosity.

	if (patch.stand.is_highlatitude_peatland_stand()) {

		if (firstTempCalc) {

			for (int i = IDX; i < IDX + NACROTELM; i++) {

				// ACROTELM

				Dz[i] = Dz_acro;
				por[i] = acro_por;
				Frac_org[i] = 0.0;
				Frac_min[i] = 0.0;	
				Frac_peat[i] = 1.0 - (acro_por + Fgas);

				Fpwp_ref[i] = peat_wp;  // Wisser et al. 2011 Very low in peatlands

				if (patch.stand.first_year == date.year && !restart) {
					Frac_water_belowpwp[i] = peat_wp;
					Frac_water[i] = acro_por - peat_wp; // i.e. initially saturated
				}

				soiltype.awc_peat[i-IDX] = (por[i] - peat_wp) * Dz[i];
			}

			for (int i = IDX + NACROTELM; i < NLAYERS; i++) {

				// CATOTELM

				Dz[i] = Dz_cato;
				por[i] = cato_por; // as in RW code
				Frac_org[i] = 0.0;
				Frac_min[i] = 0.0;
				Frac_peat[i] = 1.0 - (cato_por + Fgas);

				Fpwp_ref[i] = peat_wp;  // Wisser et al. 2011 Very low in peatlands

				if (patch.stand.first_year == date.year && !restart) {
					Frac_water_belowpwp[i] = peat_wp;
					Frac_water[i] = cato_por - peat_wp; // i.e. initially saturated
				}

				soiltype.awc_peat[i-IDX] = (por[i] - peat_wp) * Dz[i];
			}

			if (mixedl > 0) { 

				Dz[MIDX] = Frac_water[MIDX] * maxh + Frac_ice[MIDX] * maxh; // MIN value: 50mm
				Frac_org[MIDX] = 0.0;
				por[MIDX] = 0.0; // As it's filled with litter water or ice.
				Frac_min[MIDX] = 0.0;
				Frac_peat[MIDX] = 0.0;
				Frac_water[MIDX] = 0.0;

				// So, Frac_air and Frac_water are NOT defined here
			}
		} // firstTempCalc
		else {

			if (mixedl > 0) {

				// Mixed layer when soilcode == organic
				Dz[MIDX] = Frac_water[MIDX] * maxh + Frac_ice[MIDX] * maxh;
				Frac_org[MIDX] = 0.0;
			}
		} // not firstTempCalc
	}
	else {

		// NON-WETLAND SOILS

		if (firstTempCalc && patch.stand.first_year == date.year) {

			for (int ii = IDX; ii<NLAYERS; ii++) {

				if (!iforganicsoilproperties || soiltype.soilcode == 8) {

					// No updates to the standard values for LPJ-GUESS soils
					por[ii] = soiltype.porosity;
					Frac_org[ii] = soiltype.organic_frac;
					Frac_min[ii] = soiltype.mineral_frac;
					Fpwp_ref[ii] = soiltype.water_below_wp;
					Frac_water_belowpwp[ii] = soiltype.water_below_wp;
				} 
				else {			
					// Use values updated using input of soil carbon
					por[ii] = soiltype.porosity_gridcell[ii-IDX];
					Frac_org[ii] = soiltype.org_frac_gridcell[ii-IDX];
					Frac_min[ii] = soiltype.min_frac_gridcell[ii-IDX];

					double layerdepth_upper = SOILDEPTH_UPPER / NSOILLAYER_UPPER;
					double layerdepth_lower = SOILDEPTH_LOWER / NSOILLAYER_LOWER;

					if (ii < IDX+NSOILLAYER_UPPER) {
						// before organic soils: Fpwp_ref[ii] = soiltype.water_below_wp;
						Fpwp_ref[ii] = soiltype.wp[ii-IDX] / layerdepth_upper; 
					}
					else {
						Fpwp_ref[ii] = soiltype.wp[ii-IDX] / layerdepth_lower; 
					}

					Frac_water_belowpwp[ii] = Fpwp_ref[ii]; // Water content initially set to the wilting point
				}

				Dz[ii] = Dz_soil;
				Frac_water[ii] = 0.0;
				Frac_peat[ii] = 0.0;
			}
		} 
		else if (firstTempCalc && (restart || patch.stand.clone_year == date.year)) {
			// initialize after restarts
			for (int ii = IDX; ii<NLAYERS; ii++) {

				if (!iforganicsoilproperties || soiltype.soilcode == 8) {

					// No updates to the standard values for LPJ-GUESS soils
					por[ii] = soiltype.porosity;
					Frac_org[ii] = soiltype.organic_frac;
					Frac_min[ii] = soiltype.mineral_frac;
					Fpwp_ref[ii] = soiltype.water_below_wp;
				}
				else {
					// Use values updated using input of soil carbon
					por[ii] = soiltype.porosity_gridcell[ii - IDX];
					Frac_org[ii] = soiltype.org_frac_gridcell[ii - IDX];
					Frac_min[ii] = soiltype.min_frac_gridcell[ii - IDX];

					double layerdepth_upper = SOILDEPTH_UPPER / NSOILLAYER_UPPER;
					double layerdepth_lower = SOILDEPTH_LOWER / NSOILLAYER_LOWER;

					if (ii < IDX + NSOILLAYER_UPPER) {
						// before organic soils: Fpwp_ref[ii] = soiltype.water_below_wp;
						Fpwp_ref[ii] = soiltype.wp[ii - IDX] / layerdepth_upper;
					}
					else {
						Fpwp_ref[ii] = soiltype.wp[ii - IDX] / layerdepth_lower;
					}
				}
				Frac_peat[ii] = 0.0;
			}
		}// firstTempCalc
	}
 }

void Soil::snowpack_dynamics(const double &snowdepth, const int& soilsurfaceindex, int& snow_active_layers) {

	// DESCRIPTION
	// Prvate function called daily from main soil temperature routine.
	// Update the depth of each snow layer, and updates the number of active snow layers.

	// INPUT:
	// snowdepth		- the daily snow depth (mm)
	// soilsurfaceindex	- index in the Dz array

	double toplayerdepth = 50.0; // mm (as in JSBACH)
	double fulldepth = NLAYERS_SNOW * toplayerdepth;
	
	// Simple snow layer scheme. The snow layer nearest the soil surface has variable depth.
	// We add/remove toplayerdepth (e.g. 50mm) layers as the snow deepens/melts. 
	if (snowdepth > fulldepth) {
		snow_active_layers = NLAYERS_SNOW;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 2] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 3] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 4] = snowdepth-(snow_active_layers -1)*toplayerdepth;
	}
	else if (snowdepth <= fulldepth && snowdepth > (NLAYERS_SNOW - 1)*toplayerdepth) { // 200-250
		snow_active_layers = NLAYERS_SNOW - 1;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 2] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 3] = snowdepth - (snow_active_layers - 1)*toplayerdepth;
	}
	else if (snowdepth <= (NLAYERS_SNOW - 1)*toplayerdepth && snowdepth > (NLAYERS_SNOW - 2)*toplayerdepth) { // 150-200
		snow_active_layers = NLAYERS_SNOW - 2;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 2] = snowdepth - (snow_active_layers - 1)*toplayerdepth;
	}
	else if (snowdepth <= (NLAYERS_SNOW - 2)*toplayerdepth && snowdepth > (NLAYERS_SNOW - 3)*toplayerdepth) { // 100-150
		snow_active_layers = NLAYERS_SNOW - 3;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = snowdepth - (snow_active_layers - 1)*toplayerdepth;
	}
	else {
		snow_active_layers = NLAYERS_SNOW - 4; // 0-100 - usually a single layer
		Dz[soilsurfaceindex - snow_active_layers] = snowdepth;
	}

}


bool Soil::soil_temp_multilayer(const double &dailyairtemp) {

	// DESCRIPTION
	// Main soil temperature routine. Called daily from dailyaccounting_patch
	// To update the temperature of each snow, soil and padding layer

	// INPUT:
	// dailyairtemp - the daily air temperature

	int layer0, layer;

	// Layer numbers
	double surf_temp;   // surface temperature [deg. C]
	double Dsnow;		// diffusivity of snow [mm2 d-1]
	double Csnow;		// heat capacity of snow [J mm-3 K-1]
	double Ksnow;		// thermal conductivity of snow [J d-1 mm-1 K-1]

	bool analyticalSolutionTest = false;

	const double ERROR_TEMP = 80.0; // degrees C


	// INITIALISE
	// ------------------------------------------

	int daynum = date.day;

	if (firstTempCalc && patch.stand.first_year == date.year) {

		for (int i = 0; i<NLAYERS; i++) {
			Frac_ice[i] = 0.0;
			// Set the soil temperature to 0 degC initally.
			// Could also initialise T_soil with Tm (MAAT) * depth(m) * 0.015 as in Wisser et al. (2011)
			T_soil[i] = 0.0;
			T_soil_yesterday[i] = 0.0;
			T_old[i] = 0.0;
			Di[i] = 0.04 * MM2_PER_M2;
			Ci[i] = Cp_air * M3_PER_MM3;
			Ki[i] = Di[i] * Ci[i];
			Dz[i] = 100.0;
			T[i] = 0.0;
		}

		// Set up padding - use calc_padding to calculate the thickness
		// of the padding layers, initially using 100mm as the lower layer, 
		// and initialise the padding layer temperatures to zero.

		calc_padding(Dz_cato, pad_dz);

		for (layer = 0; layer<PAD_LAYERS; layer++) {
			pad_temp[layer] = 0.0;
		}

		// The log is used to speed things up. Without the log you need
		// to use an exponential, with it you just multiply.
		lKorg	= log(Korg);
		lKpeat	= log(Kpeat);
		lKmin	= log(Kmin);
		lKwater	= log(Kwater);
		lKice	= log(Kice);
		lKair	= log(Kair);

		if (patch.stand.is_highlatitude_peatland_stand())
			ngroundl = NACROTELM + NCATOTELM; // Usually 15
		else
			ngroundl = NSOILLAYER;

	} // firstTempCalc
	else if (firstTempCalc && (restart || patch.stand.clone_year == date.year)) {

		// The log is used to speed things up. Without the log you need
		// to use an exponential, with it you just multiply.
		lKorg = log(Korg);
		lKpeat = log(Kpeat);
		lKmin = log(Kmin);
		lKwater = log(Kwater);
		lKice = log(Kice);
		lKair = log(Kair);

		if (patch.stand.is_highlatitude_peatland_stand())
			ngroundl = NACROTELM + NCATOTELM; // Usually 15
		else
			ngroundl = NSOILLAYER;
	}


	// Error check
	if (!valid_layer_num(ngroundl)) {
		dprintf("Soil::soil_temp_multilayer - INVALID LAYER ERROR!!!\n");
		return false;
	}

	// Reset ice_frac
	for (int i = 0; i<NLAYERS; i++) {
		Di[i] = 0.04 * MM2_PER_M2;
		Ci[i] = Cp_air * M3_PER_MM3;
		Ki[i] = Di[i] * Ci[i];
	}

	// Allocate yesterday's temperature and ice fraction to today's 
	// ------------------------------------------

	if (!firstTempCalc || restart || patch.stand.clone_year == date.year) update_from_yesterday();

	// SNOW DENSITY and THERMAL PROPRTIES
	// ------------------------------------------

	bool snow_active_old = snow_active;
	update_snow_properties(daynum, dailyairtemp, Dsnow, Csnow, Ksnow);

	if (analyticalSolutionTest) snow_active = 0;

	// MIXED LAYER
	// On wetlands there can be a layer that is filled with water, depending on the WTP.
	// Set to 1 if this is activated
	int mixedl = 0;

	// Set indexes for first active layers:
	IDX = NLAYERS - ngroundl;							// usually 7
	MIDX = NLAYERS - ngroundl - mixedl;					// usually 7 or 6

	// SNOW LAYERS
	// ------------------------------------------

	// Determine how many snow layers are active, based on snow_active, snow depth and snow_active_layers (<=NLAYERS_SNOW);
	int snow_active_layers_old = snow_active_layers;

	double snowdepth = snowpack / (snowdens / water_density); // mm
	int soilsurfaceindex = NLAYERS - ngroundl - mixedl;

	if (ifmultilayersnow && !iftwolayersoil) {

		// Allow multiple layers?
		if (snow_active) {
			snowpack_dynamics(snowdepth, soilsurfaceindex, snow_active_layers);
		}
		else {
			snow_active_layers = 0;
		}
	}
	else {

		// ONE snow layer
		if (snow_active) {

			snow_active_layers = 1; // a single layer
			Dz[soilsurfaceindex - snow_active_layers] = snowdepth;
		}
		else {
			snow_active_layers = 0;
		}
	}

	// SNOW INDEX
	SIDX = NLAYERS - ngroundl - mixedl - snow_active_layers;

	// Update the thermal properties in the snow layers
	for (int j = SIDX; j< SIDX + snow_active_layers; j++) {
		Di[j] = Dsnow;
		// correct units
		Ci[j] = Csnow;
		Ki[j] = Ksnow;
	}

	// Set fractions of soil parameters for each soil layer & litter
	// ------------------------------------------

	if (firstTempCalc && patch.stand.first_year == date.year) {
		init_hydrology_variables(); // Initialise whc[], alwhc[], and alwhc_init[] for this patch
	}

	update_layer_fractions(daynum, mixedl, MIDX);

	// Update surface (air) temperature 
	surf_temp = dailyairtemp;

	// Update temperature in SNOW layers 
	// ------------------------------------------

	if (firstTempCalc && patch.stand.first_year == date.year) {
		SIDX_old = SIDX;
		snow_days_prev = 365;
		snow_days = 0;
		snow_active = false;
		snow_active_layers = 0;
	}

	// Counting snow days for the Wania et al. compaction scheme
	if (snow_active) {

		// If there's still snow...
		snow_days++;

		if (!snow_active_old) {
			// Reset
			snow_days = 1;
		}
	}
	else {
		if (snow_active_old)
			snow_days_prev = snow_days; // The day after all the snow has melted

		snow_days = 0;
	}

	SIDX_old = SIDX;

	// Now update snow layer temperatures if snow is active && the number of snow layers has changed
	if (snow_active && snow_active_layers != snow_active_layers_old) {

		if (snow_active_layers > snow_active_layers_old) {
			// new snow layer(s) - set new layer temperatures to air temp
			for (int j = SIDX; j<SIDX + snow_active_layers - snow_active_layers_old; j++) {
				T_soil[j] = surf_temp;
			}
			// The existing snow layers keep their temperatures
		}
		else if (snow_active_layers < snow_active_layers_old) {
			T_soil[SIDX - 1] = surf_temp;
		}
	}

	// AIR layer
	// ------------------------------------------

	Dz[SIDX - 1] = AIR_THICKNESS; // mm
	Frac_air[SIDX - 1] = 1.0;
	por[SIDX - 1] = 1.0;
	Frac_org[SIDX - 1] = 0.0;
	Frac_min[SIDX - 1] = 0.0;
	Frac_peat[SIDX - 1] = 0.0;
	Di[SIDX - 1] = Kair / Cp_air * SECS_PER_DAY * MM2_PER_M2;	// mm2 day-1
	Ci[SIDX - 1] = Cp_air * M3_PER_MM3;					// unit conversion J mm-3 K-1
	Ki[SIDX - 1] = Di[SIDX - 1] * Ci[SIDX - 1];			// units: (mm2 day-1) * (J mm-3 K-1) = J day-1 mm-1 K-1 , as required

	// Surface air temperature
	T_soil[SIDX - 1] = surf_temp;

	// NB - layer0 defined here
	layer0 = SIDX - 1;

	// Initialise unused T_soil values.
	for (int i = 0; i<SIDX - 1; i++) T_soil[i] = MISSING_VALUE;

	// Update the heat capacities, thermal conductivities and diffusivities (Di[i]) in the soil layers
	if (!iftwolayersoil)
		update_soil_diffusivities(daynum, analyticalSolutionTest);

	// Do one timestep of the Crank-Nicholson method.
	// ------------------------------------------

	if (firstTempCalc && patch.stand.first_year == date.year) calc_padding(Dz[NLAYERS - 1], pad_dz);

	// Initialise T array for the layers
	for (int i = SIDX - 1; i<NLAYERS; i++) T[i] = T_soil[i];

	if (!iftwolayersoil) {

		// Calculate the temperature - added the timesteps loop, just in case
		for (int timestep = 0; timestep < TIMESTEPS; timestep++) {

			double dayfrac = Dt / (double)TIMESTEPS;

			// Wania at al. (2009a) algorithm. Assumes vertical homogeneity in Ki and Ci.
			cnstep(layer0, Di, Dz, surf_temp, dayfrac, pad_dz, T, pad_temp);
			// Alternative algorithm that does not assume vertical homogeneity in Ki and Ci:
			// cnstep_full(layer0, Di, Dz, surf_temp, dayfrac, pad_dz, T, pad_temp, Ki, Ci);

		}
	}

	// ERROR Checking!
	for (int i = SIDX; i<NLAYERS; i++) {

		if (Di[i] <= 0.0)
			dprintf("%g%s%g\n", (double)i, ", BAD Di[i]: ", Di[i]);

		if (T[i] < -ERROR_TEMP || T[i] > ERROR_TEMP) { // -70 exceeded in Siberia sometimes

			dprintf("Soil::soil_temp_multilayer - SOIL TEMP. ERROR!!!\n");
			dprintf("%g%s%g\n", (double)date.year, "  ", (double)date.day);
			dprintf("%g%s%g%s%g\n", (double)i, ":  T[i]: ", T[i], ", Di[i]: ", Di[i]);

			return false;
		} 
		else 
			T_soil[i] = T[i]; // Store temperature values
	}

	// PHASE CHANGES
	// ------------------------------------------

	// - Global Boolean to regulate phase changes
	// - Added 100-year (approx) condition to get the water balance right during spin-up
	if (!iftwolayersoil && ifallowphasechanges && date.year >= FIRST_FREEZE_YEAR) {
		update_ice_fraction(daynum, MIDX);
	} // allowPhaseChanges?

	// Update whc[], awhc[], and wcont[]
	if (!update_layer_water_content(daynum))
		return false;

	if (analyticalSolutionTest) {
		for (int i = SIDX - 1; i < NLAYERS; i++) Frac_water[i] = 0.0; // No water if analytic soln. sought;
	}

	// THAW DEPTH
	// ------------------------------------------

	if (daynum == 0)
		maxthawdepththisyear = 0.0;

	// Define the thaw depth as the 0degC isotherm (Muller 1947, Burn 1998), 
	// which equals the depth, until which liquid water exists continuously from the top down.

	double thaw_depth_today;	// daily thawing depth, where at least some water is found, but ice can be present [mm]
	double thaw_top;	// thawing depth, where all the ice has melted [mm]

	thaw_depth_today = 0.0;

	if (!iftwolayersoil) {

		// Include padding layers
		int i = IDX;
		while (i < NLAYERS && T_soil[i] > 0.0) {
			thaw_depth_today += Dz[i];
			i++;
		}
	} 
	else {
		thaw_depth_today = SOILDEPTH_UPPER + SOILDEPTH_LOWER; // no freezing
	}

	// Record the daily thawdepth. NB! Needed for fire. 
	thaw = thaw_depth_today;

	if (thaw > maxthawdepththisyear)
		maxthawdepththisyear = thaw;

	// Calculate also at which depth there is no ice anymore. This is
	// needed in the waterbalance SR.

	thaw_top = 0.0;
	int j = IDX;
	while (j < NLAYERS && Frac_ice[j]<0.001) {
		thaw_top += Dz[j];
		j++;
	}

	if (!iftwolayersoil) {

		// Record Frac_ice and T_soil EVERY day now, not just Dec 31
		for (int i = 0; i<NLAYERS; i++) {
			Frac_ice_yesterday[i] = Frac_ice[i];
			T_soil_yesterday[i] = T_soil[i];
		}

		// Finally, update soil.temp25 with the calculated value:
		int soil25_ix = IDX + 2;
		if (Dz_soil == 100.0) {
			temp25 = T_soil[soil25_ix];
		}

		// Update monthly, later averages
		int l_ix = IDX;
		for (int sl = 0; sl < SOILTEMPOUT; sl++) {
			T_soil_monthly[date.month][sl] += T_soil[l_ix + sl] / (double)date.ndaymonth[date.month];
		}
	}

	// Not the first time this method is called	
	firstTempCalc = false;

	return true;

} // Soil::calcsoiltemp

double Soil::nmass_avail(int pref) {
	double nmass = 0.0;
	if (!ifntransform) {
		pref = NH4;
		if (NO3_mass > 0.0) {
			NH4_mass += NO3_mass;
			NO3_mass = 0.0;
		}
	}
	if (pref == NO) {
		nmass = NH4_mass + NO3_mass;
	} else if (pref == NH4) {
		nmass = NH4_mass;
	} else if (pref == NO3) {
		nmass = NO3_mass;
	}
	return nmass;
}

void Soil::nmass_subtract(double nmass, int pref) {

	if (!ifntransform) {
		pref = NH4;
	}
	double residual = 0.0;
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		if (nmass_tot > 0.0) {
			residual = nmass_tot - nmass;
			if (residual >= 0.0 || negligible(residual,10)) {
				nmass_subtract(nmass * NH4_mass / nmass_tot,NH4);
				nmass_subtract(nmass * NO3_mass / nmass_tot,NO3);
			} else {
				fail("tried to subtract more N than available");
			}
		}
	} else if (pref == NH4) {
		if (nmass<NH4_mass || negligible(nmass-NH4_mass,10)) {
			NH4_mass -= nmass;
		} else {
			double r = nmass - NH4_mass;
			NH4_mass = 0.0;
			patch.fluxes.report_flux(Fluxes::NH3_SOIL, r);
			dprintf("NH4 mass %f N substr %f \n",NH4_mass,nmass);
			//fail("tried to subtract more NH4 than available");
		}
	} else if (pref == NO3) {
		if (nmass<NO3_mass || negligible(nmass-NO3_mass,10)) {
			NO3_mass -= nmass;
		} else {
			double r = nmass - NO3_mass;
			NO3_mass = 0.0;
			patch.fluxes.report_flux(Fluxes::NO_SOIL, r);
			dprintf("NO3 mass %f N substr %f \n",NO3_mass,nmass);
			//fail("tried to subtract more NO3 than available");
		}
	}

}

void Soil::nmass_inc(double nmass, int pref) {
	if (nmass < 0.0) {
		nmass_subtract(fabs(nmass), pref);
	} else {
		if (!ifntransform) {
			pref = NH4;
		}
		if (pref == NO) {
			double nmass_tot = NH4_mass + NO3_mass;
			NH4_mass+=nmass * NH4_mass / nmass_tot;
			NO3_mass+=nmass * NO3_mass / nmass_tot;
		} else if (pref == NH4) {
			NH4_mass += nmass;
		} else if (pref == NO3) {
			NO3_mass += nmass;
		}
	}
}

void Soil::nmass_multiplic_inc(double inc, int pref) {
	if (!ifntransform) {
		pref = NH4;
	}
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		NH4_mass *= inc * NH4_mass / nmass_tot;
		NO3_mass *=inc * NO3_mass / nmass_tot;
	} else if (pref == NH4) {
		NH4_mass *= inc;
	} else if (pref == NO3) {
		NO3_mass *= inc;
	}
}

// serialize new soil variables, when finalised
void Soil::serialize(ArchiveStream& arch) {
	arch & wcont
		& wcont_evap
		& awcont_upper
		& dwcontupper
		& dwcontlower
		& mwcontupper
		& mwcontlower
		& mwcont
		& snowpack
		& snow_active
		& snow_days
		& snow_days_prev
		& snow_active_layers
		& snow_water
		& snow_ice
		& msnowdepth
		& dsnowdepth
		& thaw
		& runoff
		& temp25
		& Dz // optimise through init after restart?
		& T_old
		& T_soil_yesterday
		& T_soil_monthly
		& dtemp
		& mtemp
		& gtemp
		& pad_temp
		& pad_dz
		& cpool_slow
		& cpool_fast
		& decomp_litter_mean
		& k_soilfast_mean
		& k_soilslow_mean
		& alag
		& exp_alag
		& Frac_ice_yesterday
		& Frac_water_belowpwp
		& Frac_water
		& Frac_air
		& whc
		& alwhc
		& aw_max // optimise through init after restart?
		// Methane parameters
		// Some of these could possibly be removed to 
		// optimise for memory in state files
		& awtp
		& wtp
		& ch4_store
		& co2_store
		& CO2_soil_yesterday
		& CH4_yesterday
		& CH4_diss_yesterday
		& CH4_gas_yesterday
		& CH4_gas_vol
		& O2;
		
	for (int i = 0; i<NSOMPOOL; i++) {
		arch & sompool[i];
	}

	arch & dperc
		& orgleachfrac
		& NO2_mass
		& NO_mass
		& N2O_mass
		& N2_mass
		& NH4_mass
		& NO3_mass
		& NH4_input
		& NO3_input
		& anmin
		& animmob
		& aminleach
		& aorgNleach
		& aorgCleach
		& anfix
		& anfix_calc
		& anfix_mean
		& snowpack_NH4_mass
		& snowpack_NO3_mass
		& solvesomcent_beginyr
		& solvesomcent_endyr
		& solvesom
		& fnuptake_mean
		& morgleach_mean
		& mminleach_mean
		& labile_carbon
		& pH;
}


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// References not related to Arctic and wetland code:
//
// Fukusako, S. Int J Thermophys (1990) 11: 353. doi:10.1007/BF01133567
// Chadburn, S., Burke, E., Essery, R., Boike, J., Langer, M., Heikenfeld, M., Cox, P., 
//   and Friedlingstein, P.: An improved representation of physical permafrost dynamics 
//   in the JULES land - surface model, Geosci.Model Dev., 8, 1493 - 1508, 
//   https ://doi.org/10.5194/gmd-8-1493-2015, 2015
// Jobbagy, E. G. & Jackson, R. B. (2000) The vertical distribution of soil carbon
//   and its relation to climate and vegetation. Ecological Applications 10(2): 423-436
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Levine, J. S. (1996) Biomass Burning and Global Change. Remote Sensing, Modeling
//   and Inventory Development, and Biomass Burning in Africa, 1J. S. Levine,
//   XXXV-XLIII, MIT Press, Mass.
//
///////////////////////////////////////////////////////////////////////////////////////
//  Arctic and wetland references:
//
// Aerts, R., Verhoeven, J.T.A., & Whigham, D.F. (1999) Plant-mediated controls on nutrient cycling 
//   in temperate fenns and bogs, Ecology, 80 (7), 2170-2181. 
// Best et al. (2011) Geosci. Model Dev., 4, 677-699. www.geosci-model-dev.net/4/677/2011/
//   doi:10.5194/gmd-4-677-2011
// Bonan, G.B. 2002 Ecological Climatology, Cambridge University Press 
// Clein, J. S., and J. P. Schimel, Microbial activity of tundra and taiga soils at sub-zero 
//   temperatures, Soil Biol. Biochem., 27(9), 1231-1234, 1995.
// Cronk, J. K. and Fennessy, M. S.: Wetland Plants: Biology and Ecology, CRC Press LLC, 2001.
// Ekici, A., Chadburn, S., Chaudhary, N., Hajdu, L. H., Marmy, A., Peng, S., Boike, J., Burke, E., 
//   Friend, A. D., Hauck, C., Krinner, G., Langer, M., Miller, P. A., and Beer, C.: Site-level model 
//   intercomparison of high latitude and high altitude soil thermal dynamics in tundra and barren 
//   landscapes, The Cryosphere, 9, 1343-1361, https://doi.org/10.5194/tc-9-1343-2015, 2015.
// Frolking, S., Roulet, N. T., Moore, T. R., Richard, P. J. H., Lavoie, M., and Muller, S. D.: Modeling 
//   northern peatland decomposition and peat accumulation, Ecosystems, 4, 479-498, 2001.
// Frolking, S., Roulet, N. T., Tuittila, E., Bubier, J. L., Quillet, A., Talbot, J., and Richard, 
//   P. J. H.: A new model of Holocene peatland net primary production, decomposition, water balance, 
//   and peat accumulation, Earth Syst. Dynam., 1, 1-21, https://doi.org/10.5194/esd-1-1-2010, 2010.
// Fukusako, S. (1990) Thermophysical properties of ice, snow, and sea ice. Int. J. Thermophys., 11(2):
//	 353-372. doi:10.1007/BF01133567
// Granberg, et al. 1999 A simple model for simulation of water content, soil frost, 
//   and soil temperatures in boreal mixed mires. Water Resour. Res., 35(12), 3771-3782.
// Hillel D., (1982) Introduction to soil physics. Academic Press, San Diego, CA, USA.
// Ise, T., Dunn, A.L, Wofsy, S.C. & Moorcroft, P.R.: High sensitivity of peat decomposition to climate 
//   change through water-table feedback. Nature Geoscience volume 1, pages 763-766 (2008)
// Jaehne, B., Heinz, G., and Dietrich, W.: Measurement of the diffusion coefficients of sparingly 
//   soluble gases in water, J. Geophys. Res., 92, 10767-10776, 1987.
// Koven, C.D, Ringeval, B., Friedlingstein, P., Ciais, P., Cadule, P., Khvorostyanov, D.,
//   Krinner, G., and Tarnocai C. (2011) Permafrost carbon-climate feedbacks accelerate global warming, 
//   PNAS, vol. 108 no. 36, 14769-14774. 
// Lawrence, D. M., and A. G. Slater, 2008: Incorporating organic soil into a global climate model. 
//   Climate Dynamics, 30, 145-160, doi:10.1007/s00382-007-0278-1.
// Ling, F., and Zhang, T. (2006) Sensitivity of ground thermal regime and surface energy fluxes to
//   tundra snow density in northern Alaska. Cold Regions Science and Technology 44 (2006) 121-130
// McGuire, A. D., Christensen, T. R., Hayes, D., Heroult, A., Euskirchen, E., Kimball, J. S., Koven, C., 
//   Lafleur, P., Miller, P. A., Oechel, W., Peylin, P., Williams, M., and Yi, Y.: An assessment of 
//   the carbon balance of Arctic tundra: comparisons among observations, process models, and atmospheric 
//   inversions, Biogeosciences, 9, 3185-3204, https://doi.org/10.5194/bg-9-3185-2012, 2012.
// Potter, C. S., Davidson, E. A., and Verchot, L. V.: Estimation of global biogeochemical controls and 
//   seasonality in soil methane consumption, Chemosphere, 32, 2219-2246, 1996.
// Riera, J. L., Schindler, J. E., and Kratz, T. K.: Seasonal dynamics of carbon dioxide and methane in 
//   two clear-water lakes and two bog lakes in northern Wisconsin, USA, 
//   Can. J. Fish. Aquat. Sci., 56, 265-274, 1999.
// Sander, R.: Compilation of Henry's Law Constants for inorganic and organic species of potential 
//   importance in environmental chemistry, Tech. Rep. Version 3, MPI Mainz, Air Chemistry Department,
//   Max-Planck Institute of Chemistry, 1999.
// Schimel, J. P.: Plant transport and methane production as controls on methane flux from arctic wet 
//   meadow tundra, Biogeochem., 28, 183-200, 1995.
// Smolders, A. J. P., H. B. M. Tomassen, H. W. Pijnappel, L. P. M. Lamers, and J. G. M. Roelofs (2001), 
//   Substrate-derived CO2 is important in the development of Sphagnum spp., New Phytol., 152(2), 325- 332.
// Spahni, R., Wania, R., Neef, L., van Weele, M., Pison, I., Bousquet, P., Frankenberg, C., Foster, 
//   P. N., Joos, F., Prentice, I. C., and van Velthoven, P.: Constraining global methane emissions and 
//   uptake by ecosystems, Biogeosciences, 8, 1643-1665, https://doi.org/10.5194/bg-8-1643-2011, 2011.
// Sturm, M., Holmgren, J., Konig, M., and Morris, K. (1997) The thermal conductivity of seasonal snow. 
//   Journal of Glaciology 43 (143), 26-41.
// Swenson, S.C., Lawrence, D.M., and Lee, H. 2012. Improved Simulation of the Terrestrial Hydrological 
//   Cycle in Permafrost Regions by the Community Land Model. JAMES, 4, M08002. DOI:10.1029/2012MS000165.
// Tang, J., Miller, P. A., Persson, A., Olefeldt, D., Pilesjo, P., Heliasz, M., Jackowicz-Korczynski, 
//   M., Yang, Z., Smith, B., Callaghan, T. V., and Christensen, T. R.: Carbon budget estimation of a 
//   subarctic catchment using a dynamic ecosystem model at high spatial resolution, 
//   Biogeosciences, 12, 2791-2808, doi:10.5194/bg-12-2791-2015, 2015.
// Wania, R., Ross, I., & Prentice, I.C. (2009a) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: I. Evaluation and sensitivity of physical 
//   land surface processes. Global Biogeochemical Cycles, 23, GB3014, doi:10.1029/2008GB003412
// Wania, R., Ross, I., & Prentice, I.C. (2009b) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: II. Evaluation and sensitivity of vegetation 
//   and carbon cycle processes. Global Biogeochemical Cycles, 23, GB015, doi:10.1029/2008GB003413
// Wania, R., Ross, I., & Prentice, I.C. (2010) Implementation and evaluation of a new methane 
//   model within a dynamic global vegetation model: LPJ-WHyMe v1.3.1, Geosci. Model Dev., 3, 565-584.
// Wisser, D., Marchenko, S., Talbot, J., Treat, C., and Frolking, S. (2011) Soil temperature response 
//   to 21st century global warming: the role of and some implications for peat carbon in thawing 
//   permafrost soils in North America, Earth Syst. Dynam., 2, 121-138
// Wolf, A., Callaghan T.V., & Larson K. (2008) Future changes in vegetation and ecosystem 
//   function of the Barents Region. Climatic Change, 87:51-73 DOI 10.1007/s10584-007-9342-4
// Yurova, A., Wolf, A., Sagerfors, J., & Nilsson, M. (2007) Variations in net ecosystem 
//   exchange of carbon dioxide in a boreal mire: Modeling mechanisms linked to water table 
//   position, Journal of Geophysical Research, 112, art. no. G02025, doi:10.1029/2006JG000342.
// Zhang, W., Miller, P.A., Smith, B., Wania, R., Koenigk, T. & Doscher, R., 2013, 
//   Tundra shrubification and tree-line advance amplify arctic climate warming: results from an 
//   individual-based dynamic vegetation model. Environmental Research Letters 8: 034023.
