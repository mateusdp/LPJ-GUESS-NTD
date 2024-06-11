///////////////////////////////////////////////////////////////////////////////////////
/// \file somdynam.cpp
/// \brief Soil organic matter dynamics
///
/// \author Ben Smith (LPJ SOM dynamics, CENTURY), David Wårlind (CENTURY), Mateus Dantas (P-Cycle)
/// $Date: 2022-09-13 10:47:57 +0200 (Tue, 13 Sep 2022) $
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
#include "somdynam.h"
#include "ntransform.h"
#include "driver.h"
#include <assert.h>
#include <bitset>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

// Turnover times (in years, approximate) for litter and SOM fractions at 10 deg C with
// ample moisture (Meentemeyer 1978; Foley 1995)

static const double TAU_LITTER=2.85; // Thonicke, Sitch, pers comm, 26/11/01
static const double TAU_SOILFAST=33.0;
static const double TAU_SOILSLOW=1000.0;

static const double TILLAGE_FACTOR = 33.0 / 17.0;
// (Inverse of "tillage factor" in Chatskikh et al. 2009, Value selected by T.Pugh)


static const double FASTFRAC=0.985;
	// fraction of litter decomposition entering fast SOM pool
static const double ATMFRAC=0.7;
	// fraction of litter decomposition entering atmosphere

// Corresponds to the amount of soil available nitrogen where SOM C:N ratio reach
// their minimum (nitrogen saturation) (Parton et al 1993, Fig. 4)
// Comment: NMASS_SAT is too high when considering BNF - Zaehle
static const double NMASS_SAT = 0.002 * 0.05;
//static const double PMASS_SAT = 0.002 * 0.001;
static const double PMASS_SAT = 0.002;
// Corresponds to the nitrogen concentration in litter where SOM C:N ratio reach
// their minimum (nitrogen saturation) (Parton et al 1993, Fig. 4)
static const double NCONC_SAT = 0.02;
//static const double PCONC_SAT = 0.004;
static const double PCONC_SAT = 0.02;

//Phosphorus Constants
// rate constant for sorbed P [d-1] (Wang et al. 2007)
static const double USORB = 0.0067 / date.year_length();
//static const double USORB = 0.0067;

// rate constant for strongly sorbed P [d-1] (Wang et al. 2007)
static const double USSORB = 0.0067 / date.year_length();
//static const double USSORB = 0.0067;

// rate constant for occluded P [d-1] (Wang et al. 2007)
static const double UOCC = 1.0E-5 / date.year_length();

///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES

// Exponential decay constants for litter and SOM fractions
// Values set from turnover times (constants above) on first call to decayrates

static double k_litter10;
static double k_soilfast10;
static double k_soilslow10;

static bool firsttime=true;
	// indicates whether function decayrates has been called before


///////////////////////////////////////////////////////////////////////////////////////
// SETCONSTANTS
// Internal function (do not call directly from framework)

void setconstants() {

	// DESCRIPTION
	// Calculate exponential decay constants (annual basis) for litter and
	// SOM fractions first time function decayrates is called

	k_litter10=1.0/TAU_LITTER;
	k_soilfast10=1.0/TAU_SOILFAST;
	k_soilslow10=1.0/TAU_SOILSLOW;
	firsttime=false;
}


///////////////////////////////////////////////////////////////////////////////////////
// BALANCE LABILE AND SORBED P POOLS
// Internal function (do not call directly from framework)

void pmass_add(Soil &soil, double delta) {

	double a = -1.0;
	double b = -soil.pmass_labile - soil.soiltype.kplab - soil.soiltype.spmax + soil.pmass_sorbed + delta;
	double c = soil.pmass_sorbed * soil.pmass_labile + soil.pmass_sorbed  * soil.soiltype.kplab + delta * soil.pmass_labile + delta * soil.soiltype.kplab - soil.pmass_labile * soil.soiltype.spmax;

	double bha = std::max(0.0, std::pow(b, 2.0) - 4.0 * a * c);

	double labile_inc = (-b - std::sqrt(bha)) / 2.0 * a;
	double sorbed_inc = delta - labile_inc;

	soil.pmass_labile += labile_inc;
	soil.pmass_sorbed += sorbed_inc;

}


///////////////////////////////////////////////////////////////////////////////////////
// DECAYRATES
// Internal function (do not call directly from framework)
// used by som_dynamic_lpj()

void decayrates(double wcont,double gtemp_soil,double& k_soilfast,double& k_soilslow,
	double& fr_litter, double& fr_soilfast, double& fr_soilslow, bool tillage) {

	// DESCRIPTION
	// Calculation of fractional decay amounts for litter and fast and slow SOM
	// fractions given current soil moisture and temperature

	// INPUT PARAMETERS
	// wcont       = water content of upper soil layer (fraction of AWC)
	// gtemp_soil  = respiration temperature response incorporating damping of Q10
	//               response due to temperature acclimation (Eqn 11, Lloyd & Taylor
	//               1994)
	// tillage_fact	 = tillage factor scaling daily decay constant for fast SOM fraction (k_soilfast),
					// default values: 1.0 for no-till soils, 2.0 for conventially tilled soils.
					// Inverse of "tillage factor" in Chatskikh et al. 2009 (mean value 1 / 0.48).

	// OUTPUT PARAMETERS
	// k_soilfast  = adjusted daily decay constant for fast SOM fraction
	// k_soilslow  = adjusted daily decay constant for slow SOM fraction
	// fr_litter   = litter fraction remaining following today's decomposition
	// fr_soilfast = fast SOM fraction remaining following today's decomposition
	// fr_soilslow = slow SOM fraction remaining following today's decomposition

	double moist_response; // moisture modifier of decomposition rate

	// On first call only: set exponential decay constants

	if (firsttime) setconstants();

	// Calculate response of soil respiration rate to moisture content of upper soil layer
	// Foley 1995 Eqn 19

	moist_response=0.25+0.75*wcont;

	// Calculate litter and SOM fractions remaining following today's decomposition
	// (Sitch et al 2000 Eqn 71) adjusting exponential decay constants by moisture and
	// temperature responses and converting from annual to daily basis
	// NB: Temperature response (gtemp; Lloyd & Taylor 1994) set by framework

	k_soilfast=k_soilfast10*gtemp_soil*moist_response/(double)date.year_length();
	if (tillage) {
		k_soilfast *= TILLAGE_FACTOR; // Increased HR for crops (tillage)
	}

	k_soilslow=k_soilslow10*gtemp_soil*moist_response/(double)date.year_length();

	fr_litter=exp(-k_litter10*gtemp_soil*moist_response/(double)date.year_length());
	fr_soilfast=exp(-k_soilfast);
	fr_soilslow=exp(-k_soilslow);
}


///////////////////////////////////////////////////////////////////////////////////////
// DECAYRATES
// Should be called by framework on last day of simulation year, following call to
// som_dynamics, once annual litter production and vegetation PFT composition are close
// to their long term equilibrium (typically 500-1000 simulation years).
// NB: should be called ONCE ONLY during simulation for a particular grid cell

void equilsom_lpj(Soil& soil) {

	// DESCRIPTION
	// Analytically solves differential flux equations for fast and slow SOM pools
	// assuming annual litter inputs close to long term equilibrium

	// INPUT PARAMETER (class defined in framework header file)
	// soil = current soil status

	double nyear;
		// number of years over which decay constants and litter inputs averaged

	nyear = soil.soiltype.solvesom_end - soil.soiltype.solvesom_begin + 1;

	soil.decomp_litter_mean/=nyear;
	soil.k_soilfast_mean/=nyear;
	soil.k_soilslow_mean/=nyear;

	soil.cpool_fast=(1.0-ATMFRAC)*FASTFRAC*soil.decomp_litter_mean/
		soil.k_soilfast_mean;
	soil.cpool_slow=(1.0-ATMFRAC)*(1.0-FASTFRAC)*soil.decomp_litter_mean/
		soil.k_soilslow_mean;
}


///////////////////////////////////////////////////////////////////////////////////////
// SOM DYNAMICS
// To be called each simulation day for each modelled area or patch, following update
// of soil temperature and soil water.

void som_dynamics_lpj(Patch& patch, bool tillage) {

	// DESCRIPTION
	// Calculation of soil decomposition and transfer of C between litter and soil
	// organic matter pools.

	double k_soilfast; // adjusted daily decay constant for fast SOM fraction
	double k_soilslow; // adjusted daily decay constant for slow SOM fraction
	double fr_litter;
		// litter fraction remaining following one day's/one month's decomposition
	double fr_soilfast;
		// fast SOM fraction remaining following one day's/one month's decomposition
	double fr_soilslow;
		// slow SOM fraction remaining following one day's/one month's decomposition
	double decomp_litter; // litter decomposition today/this month (kgC/m2)
	double cflux; // accumulated C flux to atmosphere today/this month (kgC/m2)
	int p;

	// Obtain reference to Soil object
	Soil& soil=patch.soil;

	// Calculate decay constants and rates given today's soil moisture and
	// temperature

	decayrates(soil.get_soil_water_upper(), soil.gtemp, k_soilfast, k_soilslow, fr_litter,
		fr_soilfast, fr_soilslow, tillage);

	// From year soil.solvesom_begin, update running means for later solution
	// (at year soil.solvesom_end) of equilibrium SOM pool sizes

	if (date.year>=soil.soiltype.solvesom_begin) {
		soil.k_soilfast_mean+=k_soilfast;
		soil.k_soilslow_mean+=k_soilslow;
	}

	// Reduce litter and SOM pools, sum C flux to atmosphere from decomposition
	// and transfer correct proportions of litter decomposition to fast and slow
	// SOM pools

	// Reduce individual litter pools and calculate total litter decomposition
	// for today/this month

	decomp_litter=0.0;

	// Loop through PFTs

	for (p=0;p<npft;p++) {	//NB. also inactive pft's

		// For this PFT ...

		decomp_litter+=(patch.pft[p].cmass_litter_leaf+
			patch.pft[p].cmass_litter_root+
			patch.pft[p].cmass_litter_myco +
			patch.pft[p].cmass_litter_sap+
			patch.pft[p].cmass_litter_heart+
			patch.pft[p].cmass_litter_repr)*(1.0-fr_litter);

		patch.pft[p].cmass_litter_leaf*=fr_litter;
		patch.pft[p].cmass_litter_root*=fr_litter;
		patch.pft[p].cmass_litter_myco *= fr_litter;
		patch.pft[p].cmass_litter_sap*=fr_litter;
		patch.pft[p].cmass_litter_heart*=fr_litter;
		patch.pft[p].cmass_litter_repr*=fr_litter;
	}

	if (date.year>=soil.soiltype.solvesom_begin)
		soil.decomp_litter_mean+=decomp_litter;

	// Partition litter decomposition among fast and slow SOM pools
	// and flux to atmosphere

	// flux to atmosphere
	cflux=decomp_litter*ATMFRAC;

	// remaining decomposition - goes to ...
	decomp_litter-=cflux;

	// ... fast SOM pool ...
	soil.cpool_fast+=decomp_litter*FASTFRAC;

	// ... and slow SOM pool
	soil.cpool_slow+=decomp_litter*(1.0-FASTFRAC);

	// Increment C flux to atmosphere by SOM decomposition
	cflux+=soil.cpool_fast*(1.0-fr_soilfast)+soil.cpool_slow*(1.0-fr_soilslow);

	// Reduce SOM pools
	soil.cpool_fast*=fr_soilfast;
	soil.cpool_slow*=fr_soilslow;

	// Updated soil fluxes. In wetlands, some portion of cflux can be emitted as CH4, so save cflux as dcflux_soil until Soil::methane() is called.  
	if (patch.stand.landcover != PEATLAND)
		patch.fluxes.report_flux(Fluxes::SOILC, cflux);
	else 
		soil.dcflux_soil=cflux;

	// Solve SOM pool sizes at end of year given by soil.solvesom_end

	if (date.year==soil.soiltype.solvesom_end && date.islastmonth && date.islastday)
		equilsom_lpj(soil);

}

/////////////////////////////////////////////////
// CENTURY SOM DYNAMICS

/// Data type representing a selection of SOM pools
/** A selection of SOM pools is represented by a bitset,
 *  the selected pools have their corresponding bits switched on.
 */
typedef std::bitset<NSOMPOOL> SomPoolSelection;


/// Reduce decay rates to keep the daily nitrogen balance in the soil
/** Only a selected subset of the SOM pools (as specified by the caller),
 *  are considered for reducion of decay rates.
 *  Usually the first time is enough (decay rate reduction of litter).
 *  Can be for phosphorus or nitrogen
 */
void reduce_decay_rates(double decay_reduction[NSOMPOOL], double net_min_pool[NSOMPOOL], const SomPoolSelection& selected, double neg_nmass_avail) {

	int neg_min_pool[NSOMPOOL] = {0};	// Keeping track on which pools that are negative

	// Add up immobilization for considered pools
	double tot_neg_min = 0.0;
	for (int p = 0; p < NSOMPOOL; p++) {
		if (selected[p] && net_min_pool[p] < 0.0) {
			tot_neg_min += net_min_pool[p];
			neg_min_pool[p] = 1;
		}
	}

	// Calculate decay reduction
	double decay_red = 0.0;

	if (tot_neg_min < neg_nmass_avail) {
		// Enough to reduce decay rates for these pools to achieve a
		// net positive mineralization
		decay_red = 1.0 - (tot_neg_min - neg_nmass_avail) / tot_neg_min;
	}
	else {
		// Need to stop these pools from decaying to be able to get
		// a net positive mineralization
		decay_red = 1.0;
	}

	// Reduce decay rate for considered pools
	for (int p = 0; p < NSOMPOOL; p++) {
		if (selected[p]) {
			decay_reduction[p] = decay_red * neg_min_pool[p];
		}
	}
}

/// Set N:C ratios for SOM pools
/** Set N:C ratios for slow, passive, humus and soil microbial pools
 *  based on mineral nitrogen pool or litter nitrogen fraction (Parton et al 1993, Fig 4)
 */
void setntoc(Soil& soil, double fac, pooltype pool, double cton_max, double cton_min,
	double fmin, double fmax) {

	if (fac <= fmin)
		soil.sompool[pool].ntoc = 1.0 / cton_max;
	else if (fac >= fmax)
		soil.sompool[pool].ntoc = 1.0 / cton_min;
	else {
		soil.sompool[pool].ntoc = 1.0 / (cton_min + (cton_max - cton_min) *
			(fmax - fac) / (fmax - fmin));
	}
}

/// Set P:C ratios for SOM pools
/** Set P:C ratios for slow, passive, humus and soil microbial pools
*  based on mineral phosphorus pool or litter phosphorus fraction (Parton et al 1988, Fig 43)
*/
void setptoc(Soil& soil, double fac, pooltype pool, double ctop_max, double ctop_min,
	double fmin, double fmax) {

	if (fac <= fmin)
		soil.sompool[pool].ptoc = 1.0 / ctop_max;
	else if (fac >= fmax)
		soil.sompool[pool].ptoc = 1.0 / ctop_min;
	else {
		soil.sompool[pool].ptoc = 1.0 / (ctop_min + (ctop_max - ctop_min) *
			(fmax - fac) / (fmax - fmin));
	}
}


/// Temperature modifier for decomposition
/** Calculate decomposition temperature modifier (in range 0-1)
 *  [A(T_soil), Eqn A9, Comins & McMurtrie 1993; ET, Friend et al 1997; abiotic
 *  effect of soil temperature, Parton et al 1993, Fig 2)
 *
 *  \param temp_soil  Soil temperature at 25 cm depth
 */
double temperature_modifier(double temp_soil) {

	double temp_mod = temp_soil > 0.0 ? max(0.0, 0.0326 + 0.00351 * pow(temp_soil, 1.652) - pow(temp_soil / 41.748, 7.19)) : 0.0;

	// Include as an overide option when: MIN_DECOMP_TEMP < temp < 0 degC
	// This increases the respiration (from 0) between MIN_DECOMP_TEMP and 0 degC - cf Koven et al. 2011
	if (ifcarbonfreeze && temp_soil <= 0.0 && temp_soil >= MIN_DECOMP_TEMP && !iftwolayersoil) {

		double decomp_at_freezing_point = 0.0326; // temp_mod above when temp_soil = 0;
		bool linear_decrease_below_freezing = false;

		if (linear_decrease_below_freezing) {
			// Alternative: Linear approach (Koven et al. 2011)
			double slope = decomp_at_freezing_point / fabs(MIN_DECOMP_TEMP);
			temp_mod = slope * temp_soil + decomp_at_freezing_point; // i.e. a linear decrease from decomp_at_freezing_point at 0C to 0 at MIN_DECOMP_TEMP (-4C).
		} 
		else {
			// Default: Q10 relationship (Schaefer & Jafarov, 2016)
			double q10_freeze = 200.5; // i.e. average of 164 and 237 based on incubation of frozen soil samples (Mikan et al., 2002)
			temp_mod = decomp_at_freezing_point * pow(q10_freeze, temp_soil / 10.0);
		}
	}

	return temp_mod;
}


/// Water Modifier for Decomposition
/** Calculate decomposition moisture modifier (in range 0-1)
 *  Friend et al 1997, Eqn 53 (Parton et al 1993, Fig 2)
 *
 *  \param wfps  Water-filled pore space
 */
double moisture_modifier(double wfps) {

	double moist_mod = wfps < 60.0 ? exp((wfps - 60.0) * (wfps - 60.0) / -800.0) : 0.000371 * wfps * wfps - 0.0748 * wfps + 4.13;
	
	return moist_mod;
}

/// Calculates CENTURY instantaneous decay rates
/** Calculates CENTURY instantaneous decay rates given soil temperature,
 *  water content of upper soil layer
 */
void decayrates_century(Soil& soil, double temp_soil, double wcont_soil, bool tillage) {

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 2010, Figure 2)
	// plus Kirschbaum et al 2001 coarse woody debris decay
	// pools SURFSTRUCT,SOILSTRUCT,SOILMICRO,SURFHUMUS,SURFMICRO,SURFMETA,SURFFWD,SURFCWD,SOILMETA,SLOWSOM,PASSIVESOM
	const double K_MAX[] = {9.5e-3, 1.9e-2, 4.2e-2, 4.8e-4, 2.7e-2, 3.8e-2, 1.1e-2, 2.2e-3, 7.0e-2, 1.7e-3, 1.9e-6};

	// Modifier for effect of soil texture
	// Eqn 5, Parton et al 1993:

	// Modify decomposition below if this is a high-latitude peatland
	bool ispeatland = soil.patch.stand.is_highlatitude_peatland_stand();

	// Modify decomposition below if this is a wetland on mineral soils
	bool ismineralwetland = soil.patch.stand.is_true_wetland_stand();

	const double texture_mod = 1.0 - 0.75 * (soil.soiltype.clay_frac + soil.soiltype.silt_frac);
	const double texture_mod_peat = 1.0 - 0.75 * (soil.soiltype.clay_frac_peat + soil.soiltype.silt_frac_peat); // = 1

	// Calculate decomposition temperature modifier (in range 0-1)
	double temp_mod = temperature_modifier(temp_soil);

	// Calculate decomposition moisture modifier (in range 0-1)
	// Water Filled Pore Spaces (wfps) % water holding capacity at wilting point (wp) and saturation capacity (wsats)
	// is calculated with the help of Cosby et al 1984;
	// use Gerten equivalents here, but wfps COULD be made depth equivalent
	const double wfps = soil.wfps(0)*100.0;
	double moist_mod = moisture_modifier(wfps);
	double moist_mod_inundated_mineral = moisture_modifier(100); // 100% WFPS for wetlands on mineral soils, 0.36 approx

	// Combined moisture and temperature modifier
	
	// simple overrides for peatlands and mineral wetlands
	double moist_mod_saturated = 1.0; // no effect unless this is peatland

	if (ispeatland) {

		moist_mod = RMOIST;

		double cmass_total = 0.0;
		for (int p = 0; p < NSOMPOOL; p++) {
			cmass_total += soil.sompool[p].cmass;
		}

		double acrotelm_climit = 7.5; // kgC/m2 - Max C content in a 30cm-deep acrotelm - see Wania et al. (2009b)
		if (cmass_total > acrotelm_climit) { // take a weighted average of the aerobic and anaerobic moisture modifiers 
			moist_mod = (acrotelm_climit * RMOIST + (cmass_total-acrotelm_climit) * RMOIST_ANAEROBIC) / cmass_total;
			// moist_mod approaches a value of RMOIST_ANAEROBIC asymptotically as cmass_total increases
		}

		moist_mod_saturated = RMOIST_ANAEROBIC / moist_mod;
	}

	if (ismineralwetland)
		moist_mod = moist_mod_inundated_mineral;

	for (int p = 0; p < NSOMPOOL; p++) {

		// Calculate decay constant
		// (dC_I/dt / C_I; Parton et al 1993, Eqns 2-4)

		double k = K_MAX[p] * temp_mod * moist_mod;

		// Include effect of recalcitrance effect of lignin
		// Parton et al 1993 Eqn 2
		// Kirschbaum et al 2001 changed the exponential term
		// from 3 to 5.

		if (p == SURFSTRUCT || p == SOILSTRUCT || p == SURFFWD || p == SURFCWD) {
			k *= exp(-5.0 * soil.sompool[p].ligcfrac);
		}
		else if (p == SOILMICRO) {
			if (ispeatland)
				k *= texture_mod_peat;
			else
				k *= texture_mod;
		}

		// Increased HR for crops (tillage)

		if (tillage && (p == SURFMICRO || p == SURFHUMUS || p == SOILMICRO || p == SLOWSOM) && !ispeatland) {
			k *= TILLAGE_FACTOR;
		}

		// Reduced decomposition for the passive and slow pools in peatlands, as they are assumed to be in the catotelm
		if (p == PASSIVESOM || p == SLOWSOM) {
			k *= moist_mod_saturated; // ensures that a modifier of RMOIST_ANAEROBIC is used.
		}

		// Calculate fraction of carbon pool remaining after today's decomposition
		soil.sompool[p].fracremain = exp(-k);

		if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
			soil.sompool[p].mfracremain_mean[date.month] += soil.sompool[p].fracremain / date.ndaymonth[date.month];
		}
	}
}

/// Transfers specified fraction (frac) of today's decomposition
/** Transfers specified fraction (frac) of today's decomposition in donor pool type
  *  to receiver pool, transferring fraction respfrac of this to the accumulated CO2
 *  flux respsum (representing total microbial respiration today)
 *  Added P parameters
 */
void transferdecomp(Soil& soil, pooltype donor, pooltype receiver,
	double frac, double respfrac, double& respsum, double& nmin_actual, double& pmin_actual,
	double& nimmob, double& pimmob, double& net_min, double& net_pmin) {

	// decrement in donor carbon pool and nitrogen pools
	double cdec = soil.sompool[donor].cdec * frac;
	double ndec = soil.sompool[donor].ndec * frac;
	double pdec = soil.sompool[donor].pdec * frac;

	// associated nitrogen increment in receiver pool (Friend et al 1997, Eqn 49)
	double ninc = cdec * (1.0 - respfrac) * soil.sompool[receiver].ntoc;

	// associated phosphorus increment in receiver pool
	double pinc = cdec * (1.0 - respfrac) * soil.sompool[receiver].ptoc;

	// if increase in receiver nitrogen greater than decrease in donor nitrogen,
	// balance must be immobilisation from mineral nitrogen pool
	// otherwise balance is nitrogen mineralisation
	if (ninc > ndec) {
		nimmob += ninc - ndec;
		net_min += ndec - ninc;
	}
	else {
		nmin_actual += ndec - ninc;
		net_min += ndec - ninc;
	}

	// if increase in receiver phosphorus greater than decrease in donor phosphorus,
	// balance must be immobilisation from mineral phosphorus pool
	// otherwise balance is phosphorus mineralisation
	if (pinc > pdec) {
		pimmob += pinc - pdec;
		net_pmin += pdec - pinc;
	}
	else {
		pmin_actual += pdec - pinc;
		net_pmin += pdec - pinc;
	}

	// "Transfer" carbon and nitrogen to receiver
	soil.sompool[receiver].delta_cmass += cdec * (1.0 - respfrac);
	soil.sompool[receiver].delta_nmass += ninc;
	soil.sompool[receiver].delta_pmass += pinc;

	// Transfer microbial respiration
	respsum += cdec * respfrac;
}

/// Fluxes between the CENTURY pools, and CO2 release to the atmosphere
/** Daily or monthly fluxes between the ten CENTURY pools, and CO2 release to the atmosphere
 *  Parton et al 1993, Fig 1; Comins & McMurtrie 1993, Appendix A
 *
 *  \param ifequilsom Whether the function is called during calculation om SOM pool equilibrium,
 *                    \see equilsom(). During this stage, somfluxes shouldn't calculate decayrates
 *                    itself, and shouldn't produce output like fluxes etc.
 */
void somfluxes(Patch& patch, bool ifequilsom, bool tillage) {


	double respsum ;
	double leachsum_cmass, leachsum_nmass, leachsum_pmass;
	double nmin_actual;	// actual (not net) nitrogen mineralisation
	double pmin_actual;	// actual (not net) phosphorus mineralisation
	double nimmob;		// nitrogen immobilisation
	double pimmob;		// phosphorus immobilisation

	const double EPS = 1.0e-16;

	Soil& soil = patch.soil;

	// mineral nitrogen mass available

	const double nmin_mass = soil.nmass_avail(NH4);// + soil.NO3_mass;
	// mineral phosphorus mass available
	const double pmin_mass = soil.pmass_labile;
	
	if (date.day == 0) {
		soil.anmin = 0.0;
		soil.animmob = 0.0;
		soil.apmin = 0.0;
		soil.apimmob = 0.0;
	}

	//// Warning if soil available nitrogen is negative (if happens once or so no problem, but if it propagates through time then it is)
	////double dummy;
	//if (ifnlim) {
	//	assert(soil.NH4_mass > -EPS);
	//	//if(soil.NH4_mass < -EPS)
	//	//dummy = 1;
	//}

	//// Warning if soil available phosphorus is negative (if happens once or so no problem, but if it propagates through time then it is)
	//if (ifplim) {
	//	assert(soil.pmass_labile > -EPS);
	//}

	// Set N:C ratios for humus, soil microbial, passive and slow pool based on estimated mineral nitrogen pool
	// (Parton et al 1993, Fig 4)
 
	// ForCent (Parton 2010) values
	setntoc(soil, nmin_mass, SLOWSOM, 30.0, 15.0, 0.0, NMASS_SAT);

	setntoc(soil, nmin_mass, SOILMICRO, 15.0, 6.0, 0.0, NMASS_SAT);

	setntoc(soil, nmin_mass, SURFHUMUS, 30.0, 15.0, 0.0, NMASS_SAT);

	// Set P:C ratios for humus, soil microbial, passive and slow pool based on estimated mineral nitrogen pool
	// (Parton et al 1988, Fig 2)

	setptoc(soil, pmin_mass, SLOWSOM, 200.0, 90.0, 0.0, PMASS_SAT);

	setptoc(soil, pmin_mass, PASSIVESOM, 200.0, 20.0, 0.0, PMASS_SAT);

	setptoc(soil, pmin_mass, SOILMICRO, 80.0, 30.0, 0.0, PMASS_SAT); //Check this

    //setptoc(soil, pmin_mass, SURFHUMUS, 200.0, 90.0, 0.0, PMASS_SAT);


	if (!ifequilsom) {

		// Calculate potential fraction remaining following decay today for all pools
		// (assumes no nitrogen limitation)
		decayrates_century(soil, soil.get_soil_temp_25(), soil.get_soil_water_upper(), tillage);

	}

	// Calculate decomposition in all pools assuming these decay rates

	// Save delta carbon and nitrogen mass

	bool net_mineralization = false;
	bool net_pmineralization = false;
	int times = 0;

	int ptimes = 0;
	double decay_reduction_np[NSOMPOOL] = { 0.0 };
	double decay_reduction_n[NSOMPOOL] = { 0.0 };
	double decay_reduction_p[NSOMPOOL] = { 0.0 };
	double init_negative_nmass, init_ntoc_reduction;
	double ntoc_reduction = 0.8;
	double init_negative_pmass, init_ptoc_reduction;
	double ptoc_reduction = 0.8;

	// If necessary, the decay rates in the pools will be reduced in groups, one group
	// is reduced after each iteration in the loop below. The groups are defined by
	// how the pools feed into each other.
	SomPoolSelection reduction_groups[4];
	reduction_groups[0].set(SURFSTRUCT).set(SURFMETA).set(SURFFWD).set(SURFCWD).set(SOILSTRUCT).set(SOILMETA);
	reduction_groups[1].set(SURFMICRO);
	reduction_groups[2].set(SURFHUMUS);
	reduction_groups[3].set(SOILMICRO).set(SLOWSOM).set(PASSIVESOM);

	// If mineralization together with soil available nitrogen is negative then decay rates are decreased
	// The SOM system have five try to get a positive result, after that all pools decay rate has been
	// affected by nitrogen limitation
	while ((!net_mineralization || !net_pmineralization) && (times < 5 || ptimes < 5)) {
	//while ((!net_mineralization && !net_pmineralization) && (times < 5 && ptimes < 5)) {

		respsum = 0.0;
		nmin_actual = 0.0;
		pmin_actual = 0.0;
		nimmob = 0.0;
		pimmob = 0.0;
		leachsum_cmass = 0.0;
		leachsum_nmass = 0.0;
		leachsum_pmass = 0.0;

		// Calculate decomposition in all pools assuming these decay rates
		for (int p = 0; p < NSOMPOOL; p++) {
			
			if(ifplim)
				decay_reduction_np[p] = max(decay_reduction_n[p], decay_reduction_p[p]);
			else
				//Mateus: Choose larger decay reduction between N and P
				decay_reduction_np[p] = decay_reduction_n[p];

			soil.sompool[p].cdec = soil.sompool[p].cmass * (1.0 - soil.sompool[p].fracremain) * (1.0 - decay_reduction_np[p]);
			soil.sompool[p].ndec = soil.sompool[p].nmass * (1.0 - soil.sompool[p].fracremain) * (1.0 - decay_reduction_np[p]);
			soil.sompool[p].pdec = soil.sompool[p].pmass * (1.0 - soil.sompool[p].fracremain) * (1.0 - decay_reduction_np[p]);

			soil.sompool[p].delta_cmass = 0.0;
			soil.sompool[p].delta_nmass = 0.0;
			soil.sompool[p].delta_pmass = 0.0;
			soil.sompool[p].delta_cmass -= soil.sompool[p].cdec;
			soil.sompool[p].delta_nmass -= soil.sompool[p].ndec;
			soil.sompool[p].delta_pmass -= soil.sompool[p].pdec;
		}

		double net_min[NSOMPOOL] = {0};
		double net_pmin[NSOMPOOL] = {0};

		// Partition potential decomposition among receiver pools

		// Donor pool SURFACE STRUCTURAL

		transferdecomp(soil, SURFSTRUCT, SURFMICRO, 1.0 - soil.sompool[SURFSTRUCT].ligcfrac,
			0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFSTRUCT], net_pmin[SURFSTRUCT]);

		transferdecomp(soil, SURFSTRUCT, SURFHUMUS, soil.sompool[SURFSTRUCT].ligcfrac, 0.3,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFSTRUCT], net_pmin[SURFSTRUCT]);

		// Donor pool SURFACE METABOLIC

		transferdecomp(soil, SURFMETA, SURFMICRO, 1.0, 0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFMETA], net_pmin[SURFMETA]);

		// Donor pool SOIL STRUCTURAL

		transferdecomp(soil, SOILSTRUCT, SOILMICRO, 1.0 - soil.sompool[SOILSTRUCT].ligcfrac,
			0.55, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILSTRUCT], net_pmin[SOILSTRUCT]);

		transferdecomp(soil, SOILSTRUCT, SLOWSOM, soil.sompool[SOILSTRUCT].ligcfrac, 0.3,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILSTRUCT], net_pmin[SOILSTRUCT]);

		// Donor pool SOIL METABOLIC

		transferdecomp(soil, SOILMETA, SOILMICRO, 1.0, 0.55, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILMETA], net_pmin[SOILMETA]);

		// Donor pool SURFACE FINE WOODY DEBRIS

		transferdecomp(soil, SURFFWD, SURFMICRO, 1.0 - soil.sompool[SURFFWD].ligcfrac,
			0.76, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFFWD], net_pmin[SURFFWD]);

		transferdecomp(soil, SURFFWD, SURFHUMUS, soil.sompool[SURFFWD].ligcfrac, 0.4,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFFWD], net_pmin[SURFFWD]);

		// Donor pool SURFACE COARSE WOODY DEBRIS

		transferdecomp(soil, SURFCWD, SURFMICRO, 1.0 - soil.sompool[SURFCWD].ligcfrac,
			0.9, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFCWD], net_pmin[SURFCWD]);

		transferdecomp(soil, SURFCWD, SURFHUMUS, soil.sompool[SURFCWD].ligcfrac, 0.5,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFCWD], net_pmin[SURFCWD]);

		// Donor pool SURFACE MICROBE

		transferdecomp(soil, SURFMICRO, SURFHUMUS, 1.0, 0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFMICRO], net_pmin[SURFMICRO]);

		// Donor pool SURFACE HUMUS

		transferdecomp(soil, SURFHUMUS, SLOWSOM, 1.0, 0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFHUMUS], net_pmin[SURFHUMUS]);

		// Donor pool SLOW SOM

		// First work out partitioning coefficients (Fig 1, Parton et al 1993)
		double csp = max(0.0, 0.003 - 0.009 * soil.get_clayfrac());
		double respfrac = 0.55;
		double csa = 1.0 - csp - respfrac;

		transferdecomp(soil, SLOWSOM, SOILMICRO, csa, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SLOWSOM], net_pmin[SLOWSOM]);

		transferdecomp(soil, SLOWSOM, PASSIVESOM, csp, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SLOWSOM], net_pmin[SLOWSOM]);

		// Account for respiration flux
		// Nitrogen associated with this respiration is mineralised (Parton et al 1993, p 791)
		respsum += respfrac * soil.sompool[SLOWSOM].cdec;

		if (!negligible(soil.sompool[SLOWSOM].cmass)) {
			nmin_actual += respfrac * soil.sompool[SLOWSOM].cdec * soil.sompool[SLOWSOM].nmass / soil.sompool[SLOWSOM].cmass;
			pmin_actual += respfrac * soil.sompool[SLOWSOM].cdec * soil.sompool[SLOWSOM].pmass / soil.sompool[SLOWSOM].cmass;
		}

		// Donor pool SOIL MICROBE

		// Fraction lost to  microbial respiration (F_t, Parton et al 1993 Eqn 7)
		respfrac = max(0.0, 0.85 - 0.68 * (soil.get_clayfrac() + soil.get_siltfrac()));

		// Fraction entering passive SOM pool (Parton et al 1993, Eqn 9)
		double cap = 0.003 + 0.032 * soil.get_clayfrac();

		transferdecomp(soil, SOILMICRO, PASSIVESOM, cap, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILMICRO], net_pmin[SOILMICRO]);

		// Fraction entering slow SOM pool
		csp = 1.0 - respfrac - soil.orgleachfrac - cap;

		transferdecomp(soil, SOILMICRO, SLOWSOM, csp, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILMICRO], net_pmin[SOILMICRO]);

		// Account for respiration flux
		// nitrogen associated with this respiration is mineralised (Parton et al 1993, p 791)
		respsum += respfrac * soil.sompool[SOILMICRO].cdec;

		// Account for organic carbon leaching loss
		leachsum_cmass = soil.orgleachfrac * soil.sompool[SOILMICRO].cdec;

		if (!negligible(soil.sompool[SOILMICRO].cmass)) {
			nmin_actual += respfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].nmass / soil.sompool[SOILMICRO].cmass;
			pmin_actual += respfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].pmass / soil.sompool[SOILMICRO].cmass;

			// Account for organic nitrogen leaching loss
			leachsum_nmass = soil.orgleachfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].nmass / soil.sompool[SOILMICRO].cmass;
			// Account for organic phosphorus leaching loss
			leachsum_pmass = soil.orgleachfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].pmass / soil.sompool[SOILMICRO].cmass;
		}

		// Donor pool PASSIVE SOM

		transferdecomp(soil, PASSIVESOM, SOILMICRO, 1.0, 0.55, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[PASSIVESOM], net_pmin[PASSIVESOM]);

		// Total net mineralization
		double tot_net_min = nmin_actual - nimmob;

		// Estimate daily soil mineral nitrogen pool after decomposition
		// (negative value = immobilisation)
		if ((tot_net_min + nmin_mass + EPS >= 0.0) || !ifnlim) {

			net_mineralization = true;
		}
		else if (!ifnlim) {

			// Not minding immobilisation higher than nmass_avail during free nitrogen years
			if (date.year > freenyears) {

				// Immobilization larger than soil available nitrogen -> reduce targeted N concentration in SOM pool with flexible N:C ratios
				if (times == 0) {
					// initial reduction
					init_negative_nmass = tot_net_min + nmin_mass;
					init_ntoc_reduction = ntoc_reduction;
				}
				else {
					// trying to match needed N:C reduction
					ntoc_reduction = min(init_ntoc_reduction, pow(init_ntoc_reduction, 1.0 / (1.0 - (tot_net_min + nmin_mass) / init_negative_nmass) + 1.0));
				}

				soil.sompool[SLOWSOM].ntoc *= ntoc_reduction;
				soil.sompool[SOILMICRO].ntoc *= ntoc_reduction;
				soil.sompool[SURFHUMUS].ntoc *= ntoc_reduction;

				net_mineralization = false;
			}
			else {
				net_mineralization = true;
			}
		}
		else {

			// Immobilization larger than soil available nitrogen -> reduce decay rates
			if (times < 4) {
				reduce_decay_rates(decay_reduction_n, net_min, reduction_groups[times], tot_net_min + nmin_mass);
			}
			net_mineralization = false;
		}
		times++;
	
		
		// Total net P mineralization
		double tot_net_pmin = pmin_actual - pimmob;

		//subtract pflux up to here to calculate total pmineralization

		//Phosphorus reduce decay rates
		// (negative value = immobilisation)
		if (tot_net_pmin + pmin_mass + EPS >= 0.0) {

			net_pmineralization = true;
		}
		else if (!ifplim) {

			// Not minding immobilisation higher than nmass_avail during free nitrogen years
			if (date.year > freenyears) {

				// Immobilization larger than soil available nitrogen -> reduce targeted N concentration in SOM pool with flexible N:C ratios
				if (ptimes == 0) {
					// initial reduction
					init_negative_pmass = tot_net_pmin + pmin_mass;
					init_ptoc_reduction = ptoc_reduction;
				}
				else {
					// trying to match needed P:C reduction
					ptoc_reduction = min(init_ptoc_reduction, pow(init_ptoc_reduction, 1.0 / (1.0 - (tot_net_pmin + pmin_mass) / init_negative_pmass) + 1.0));
				}

				soil.sompool[SLOWSOM].ptoc *= ptoc_reduction;
				soil.sompool[SOILMICRO].ptoc *= ptoc_reduction;
				soil.sompool[SURFHUMUS].ptoc *= ptoc_reduction;

				net_pmineralization = false;
			}
			else {
				net_pmineralization = true;
			}
		}
		else {

			// Immobilization larger than soil available nitrogen -> reduce decay rates
			if (ptimes < 4) {
				reduce_decay_rates(decay_reduction_p, net_pmin, reduction_groups[ptimes], tot_net_pmin + pmin_mass);
			}
			net_pmineralization = false;
		}

		ptimes++;
 
	}


	// Most of organic leaching is retained in the ecosystem (82%, Wilcke), and mineralized (Parton 1988)
	double leachsum_pmass_retained = leachsum_pmass * 1.0;

	double leachsum_pmass_lost = leachsum_pmass * 0.0;

	// Update pool sizes

	for (int p = 0; p < NSOMPOOL; p++) {
		soil.sompool[p].cmass += soil.sompool[p].delta_cmass;
		soil.sompool[p].nmass += soil.sompool[p].delta_nmass;
		soil.sompool[p].pmass += soil.sompool[p].delta_pmass;
	}

	if (!ifequilsom) {

		// Updated soil fluxes. In wetlands, some portion of cflux can be emitted as CH4, so save cflux as dcflux_soil until Soil::methane() is called.  
		if (patch.stand.landcover != PEATLAND)
			patch.fluxes.report_flux(Fluxes::SOILC, respsum);
		else 
			soil.dcflux_soil=respsum; 

		// Sum annual organic nitrogen leaching

		soil.aorgNleach += leachsum_nmass;

		// Sum annual organic phosphorus leaching

		//soil.aorgPleach += leachsum_pmass;
		// Leached organic phosphorus is assumed to be mineralized (Parton 1988, pg. 117)
		soil.aorgPleach += leachsum_pmass_lost;

		// Sum annual organic carbon leaching

		soil.aorgCleach += leachsum_cmass;

		// Sum annuals
		soil.anmin += nmin_actual;
		soil.animmob += nimmob;
		soil.apmin += pmin_actual;
		soil.apimmob += pimmob;
	}

	// Fraction of microbial resp. is assumed to produce labile carbon
	soil.labile_carbon = respsum * frac_labile_carbon;

	// Adding mineral nitrogen to soil available pool
	double nmin_inc = nmin_actual - nimmob; 
		
	//double pmin_inc = pmin_actual - pimmob;
	// Leached organic phosphorus is assumed to be mineralized (Parton 1988, pg. 117)
	double pmin_inc = pmin_actual - pimmob + leachsum_pmass_retained;

	//soil.pmass_labile = max(0.0, soil.pmass_labile + pmin_inc);
	//soil.pmass_labile_delta += pmin_inc;

	pmass_add(soil, pmin_inc);


	///////////////////////////// Balanced dynamics of P labile and P sorbed, also flux into strongly sorbed pool

	double delta_strongly_sorbed = USORB * soil.pmass_sorbed - USSORB * soil.pmass_strongly_sorbed;

	pmass_add(soil, -delta_strongly_sorbed);

	patch.fluxes.report_flux(Fluxes::P_SOIL, delta_strongly_sorbed);


	////////////////////////////////////////////////////////////////

	// Estimate of N flux from soil (simple CLM-CN approach)
	if(!ifntransform) {
		double nflux = nmin_inc > 0.0 ? (nmin_inc) * 0.01 : 0.0;

		nmin_inc -= nflux;

		if (!ifequilsom) {
			patch.fluxes.report_flux(Fluxes::NH3_SOIL, nflux);
		}
	}

	soil.nmass_inc(nmin_inc,NH4);

	// If no nitrogen limitation or during free nitrogen years set soil
	// available nitrogen to its saturation level.
	if (!ifnlim || date.year <= freenyears) {
		if(ifntransform) {
			soil.NH4_mass = NMASS_SAT / 2.0;
			soil.NO3_mass = NMASS_SAT / 2.0;
		} 
		else {
			soil.NH4_mass = NMASS_SAT;
		}
	}

	// If no phosphorus limitation or during free phosphorus years set soil
	// available phosphorus to its saturation level.
	if (!ifplim || date.year <= freenyears) {
		soil.pmass_labile = PMASS_SAT;
		//soil.pmass_labile = 0.0;
		//Value of sorbed phosphorus pool, based on labile p and soil parameters [KgP/m-2].
		//Equilibrated instantanously based on Wang 2007, 2010
		//soil.pmass_sorbed = (PMASS_SAT * soil.soiltype.spmax) / (soil.soiltype.kplab + PMASS_SAT);
		//soil.pmass_sorbed = (PMASS_SAT * soil.soiltype.spmax * soil.soiltype.kplab) / pow(soil.soiltype.kplab + PMASS_SAT, 2.0);
		soil.pmass_sorbed = soil.soiltype.spmax;
		//soil.pmass_sorbed = 0.0;
	}
}

/// Litter lignin to N ratio (for leaf and root litter)
/** Specific lignin fractions for leaf and root
 *  are specified in transfer_litter()
 *  Check if needed function for P too
 */
double lignin_to_n_ratio(double cmass_litter, double nmass_litter, double LIGCFRAC, double cton_avr) {

	if (!negligible(nmass_litter) && ifnlim) {
		return max(0.0, LIGCFRAC * cmass_litter / nmass_litter);
	}
	else {
		return max(0.0, LIGCFRAC * cton_avr / (1.0 - nrelocfrac));
	}
}

/// Metabolic litter fraction (for leaf and root litter)
/** Fm, Parton et al 1993, Eqn 1:
 *  NB: incorrect/out-of-date intercept and slope given in Eqn 1; values used in
 *  code of CENTURY 4.0 used instead (also correct in Parton et al. 1993, figure 1)
 *
 *  Check if needed function for P too
 * \param lton  Litter lignin:N ratio
 */
double metabolic_litter_fraction(double lton) {
	return max(0.0, 0.85 - lton * 0.013);
}

/// Transfers litter from growth, mortality and fire
/** Called daily to transfer last year's litter from vegetation
 *  (turnover, mortality and fire) to soil litter pools.
 *  Alternatively, with daily allocation and harvest/turnover,
 *  the litter produced a certain day.
 */
void transfer_litter(Patch& patch) {

	Soil& soil = patch.soil;

	double lat = patch.get_climate().lat;

	double EPS = -1.0e-16;

	double ligcmass_new, ligcmass_old;

	// Fire (GlobFIRM)
	double litterme[NSOMPOOL]  = {0.};
	double fireresist[NSOMPOOL]= {0.};
	if ( firemodel == GLOBFIRM ) {

		litterme[SURFSTRUCT]   = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].litterme;
		litterme[SURFMETA]     = soil.sompool[SURFMETA].cmass   * soil.sompool[SURFMETA].litterme;
		litterme[SURFFWD]      = soil.sompool[SURFFWD].cmass    * soil.sompool[SURFFWD].litterme;
		litterme[SURFCWD]      = soil.sompool[SURFCWD].cmass    * soil.sompool[SURFCWD].litterme;

		fireresist[SURFSTRUCT] = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].fireresist;
		fireresist[SURFMETA]   = soil.sompool[SURFMETA].cmass   * soil.sompool[SURFMETA].fireresist;
		fireresist[SURFFWD]    = soil.sompool[SURFFWD].cmass    * soil.sompool[SURFFWD].fireresist;
		fireresist[SURFCWD]    = soil.sompool[SURFCWD].cmass    * soil.sompool[SURFCWD].fireresist;
	}

	double leaf_litter = 0.0;
	double root_litter = 0.0;
	double myco_litter = 0.0;
	double wood_litter = 0.0;

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& pft=patch.pft.getobj();

		// For stands with yearly growth, drop leaf and root litter on first month of the year for northern hemisphere
		// and first month of the second half of the year for southern hemisphere for summergreen trees. For evergreens
		// as a fraction every day and for raingreens on the drierst month from last year. 
		// For stands with daily growth, harvest and/or turnover, do this when patch.is_litter_day is true.

		double cmass_litter_leaf, nmass_litter_leaf, pmass_litter_leaf;
		double cmass_litter_root, nmass_litter_root, pmass_litter_root;
		double cmass_litter_myco;

		double frac_lr = 0.0;
		//  Is a litter day, drop all leaf litter (crop)
		if (patch.is_litter_day) {
			frac_lr = 1.0;
		}
		// For summergreens drop leaf litter over all days during Jan in NH and July SH
		// and raingreens on the month with lowest phen
		else if ((pft.pft.phenology == SUMMERGREEN && ((lat >= 0.0 && date.month == 0) || (lat < 0.0 && date.month == 6))) ||
			(pft.pft.phenology == RAINGREEN && date.month == pft.driest_mth)) {
			frac_lr = 1.0 / (date.ndaymonth[date.month] - date.dayofmonth);
		}
		// Evergreens drops leaf litter every day
		else if (pft.pft.phenology == EVERGREEN || pft.pft.phenology == ANY) {
			frac_lr = 1.0 / (date.year_length() - date.day);
		}

		cmass_litter_leaf = pft.cmass_litter_leaf * frac_lr;
		nmass_litter_leaf = pft.nmass_litter_leaf * frac_lr;
		pmass_litter_leaf = pft.pmass_litter_leaf * frac_lr;
		cmass_litter_root = pft.cmass_litter_root * frac_lr;
		nmass_litter_root = pft.nmass_litter_root * frac_lr;
		pmass_litter_root = pft.pmass_litter_root * frac_lr;
		cmass_litter_myco = pft.cmass_litter_myco * frac_lr;

		pft.cmass_litter_leaf -= cmass_litter_leaf;
		pft.nmass_litter_leaf -= nmass_litter_leaf;
		pft.pmass_litter_leaf -= pmass_litter_leaf;
		pft.cmass_litter_root -= cmass_litter_root;
		pft.nmass_litter_root -= nmass_litter_root;
		pft.pmass_litter_root -= pmass_litter_root;
		pft.cmass_litter_myco -= cmass_litter_myco;

		// LEAF

		if (!negligible(cmass_litter_leaf) || !negligible(nmass_litter_leaf) || !negligible(pmass_litter_leaf)) {
		//if (!negligible(cmass_litter_leaf) || !negligible(nmass_litter_leaf)) {

			// Calculate inputs to surface structural and metabolic litter

			// Leaf litter lignin:N ratio
			double leaf_lton;

			if (!ifslavary || negligible(pft.cmass_leaf) || negligible(pft.nmass_leaf))
				leaf_lton = lignin_to_n_ratio(cmass_litter_leaf, nmass_litter_leaf, LIGCFRAC_LEAF, pft.pft.cton_leaf_avr); //maybe here is pft.pft.cmass_litter_leaf/pft.pft.nmass_litter_leaf
			else
				//leaf_lton = lignin_to_n_ratio(cmass_litter_leaf, nmass_litter_leaf, LIGCFRAC_LEAF, pft.cmass_litter_leaf / pft.nmass_litter_leaf);
				leaf_lton = lignin_to_n_ratio(cmass_litter_leaf, nmass_litter_leaf, LIGCFRAC_LEAF, pft.cmass_leaf / pft.nmass_leaf);

			// Metabolic litter fraction for leaf litter
			double fm = metabolic_litter_fraction(leaf_lton);

			ligcmass_old = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].ligcfrac;

			// Add to pools
			soil.sompool[SURFSTRUCT].cmass += cmass_litter_leaf * (1.0 - fm);
			soil.sompool[SURFSTRUCT].nmass += nmass_litter_leaf * (1.0 - fm);
			soil.sompool[SURFSTRUCT].pmass += pmass_litter_leaf * (1.0 - fm);
			soil.sompool[SURFMETA].cmass += cmass_litter_leaf * fm;
			soil.sompool[SURFMETA].nmass += nmass_litter_leaf * fm;
			soil.sompool[SURFMETA].pmass += pmass_litter_leaf * fm;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter(cmass_litter_leaf * (1.0 - fm), nmass_litter_leaf * (1.0 - fm), pmass_litter_leaf * (1.0 - fm), SURFSTRUCT);
				soil.litterSolveSOM.add_litter(cmass_litter_leaf * fm, nmass_litter_leaf * fm, pmass_litter_leaf * fm, SURFMETA);

			}

			// Fire
			if (firemodel == GLOBFIRM) {

				litterme[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * pft.pft.litterme;
				fireresist[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * pft.pft.fireresist;

				litterme[SURFMETA] += cmass_litter_leaf * fm * pft.pft.litterme;
				fireresist[SURFMETA] += cmass_litter_leaf * fm * pft.pft.fireresist;
			}
			// NB: reproduction litter cannot contain nitrogen!!

			ligcmass_new = cmass_litter_leaf * (1.0 - fm) * LIGCFRAC_LEAF;

			if (negligible(soil.sompool[SURFSTRUCT].cmass)) {
				soil.sompool[SURFSTRUCT].ligcfrac = 0.0;
			}
			else {
				soil.sompool[SURFSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SURFSTRUCT].cmass;
			}
		}

		pft.cmass_leaf = 0.0;
		pft.nmass_leaf = 0.0;
		pft.pmass_leaf = 0.0;

		// ROOT (Mycorrhiza c litter goes in root litter since it has no N or P mass)

		if (!negligible(cmass_litter_root) || !negligible(nmass_litter_root) || !negligible(pmass_litter_root)) {
		//if (!negligible(cmass_litter_root) || !negligible(nmass_litter_root)) {

			// Calculate inputs to soil structural and metabolic litter

			// Root litter lignin:N ratio
			double root_lton;

			if(!ifslavary || negligible(pft.cmass_root) || negligible(pft.nmass_root))
				root_lton = lignin_to_n_ratio(cmass_litter_root + cmass_litter_myco, nmass_litter_root, LIGCFRAC_ROOT, pft.pft.cton_root_avr); //maybe here is pft.pft.cmass_litter_root/pft.pft.nmass_litter_root
			else
				//root_lton = lignin_to_n_ratio(cmass_litter_root + cmass_litter_myco, nmass_litter_root, LIGCFRAC_ROOT, pft.cmass_litter_root / pft.nmass_litter_root);
				root_lton = lignin_to_n_ratio(cmass_litter_root + cmass_litter_myco, nmass_litter_root, LIGCFRAC_ROOT, pft.cmass_root / pft.nmass_root);

			// Metabolic litter fraction for root litter
			double fm = metabolic_litter_fraction(root_lton);

			ligcmass_new = (cmass_litter_root + cmass_litter_myco) * (1.0 - fm) * LIGCFRAC_ROOT;
			ligcmass_old = soil.sompool[SOILSTRUCT].cmass * soil.sompool[SOILSTRUCT].ligcfrac;

			// Add to pools and update lignin fraction in structural pool
			soil.sompool[SOILSTRUCT].cmass += (cmass_litter_root + cmass_litter_myco) * (1.0 - fm);
			soil.sompool[SOILSTRUCT].nmass += nmass_litter_root * (1.0 - fm);
			soil.sompool[SOILSTRUCT].pmass += pmass_litter_root * (1.0 - fm);
			soil.sompool[SOILMETA].cmass += (cmass_litter_root + cmass_litter_myco) * fm;
			soil.sompool[SOILMETA].nmass += nmass_litter_root * fm;
			soil.sompool[SOILMETA].pmass += pmass_litter_root * fm;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter((cmass_litter_root + cmass_litter_myco) * (1.0 - fm), nmass_litter_root * (1.0 - fm), pmass_litter_root * (1.0 - fm), SOILSTRUCT);
				soil.litterSolveSOM.add_litter((cmass_litter_root + cmass_litter_myco) * fm, nmass_litter_root * fm, pmass_litter_root * fm, SOILMETA);

			}

			if (negligible(soil.sompool[SOILSTRUCT].cmass)) {
				soil.sompool[SOILSTRUCT].ligcfrac = 0.0;
			}
			else {
				soil.sompool[SOILSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SOILSTRUCT].cmass;
			}
		}

		pft.cmass_root = 0.0;
		pft.nmass_root = 0.0;
		pft.pmass_root = 0.0;

		// WOOD

		// Woody debris enters two woody litter pools as described in
		// Kirschbaum and Paul (2002).

		// Woody litter is dropped as a fraction every day
		double frac = 1.0 / (date.year_length() - date.day);
		double cmass_litter_sap = pft.cmass_litter_sap * frac;
		double nmass_litter_sap = pft.nmass_litter_sap * frac;
		double pmass_litter_sap = pft.pmass_litter_sap * frac;
		double cmass_litter_heart = pft.cmass_litter_heart * frac;
		double nmass_litter_heart = pft.nmass_litter_heart * frac;
		double pmass_litter_heart = pft.pmass_litter_heart * frac;

		pft.cmass_litter_sap -= cmass_litter_sap;
		pft.nmass_litter_sap -= nmass_litter_sap;
		pft.pmass_litter_sap -= pmass_litter_sap;
		pft.cmass_litter_heart -= cmass_litter_heart;
		pft.nmass_litter_heart -= nmass_litter_heart;
		pft.pmass_litter_heart -= pmass_litter_heart;

		// SAP WOOD

		if (!negligible(cmass_litter_sap) || !negligible(nmass_litter_sap) || !negligible(pmass_litter_sap)) {
		//if (!negligible(cmass_litter_sap) || !negligible(nmass_litter_sap)) {

			// Fine woody debris

			ligcmass_new = cmass_litter_sap * LIGCFRAC_WOOD;
			ligcmass_old = soil.sompool[SURFFWD].cmass * soil.sompool[SURFFWD].ligcfrac;


			// Add to fine woody pool and update lignin fraction in pool
			soil.sompool[SURFFWD].cmass += cmass_litter_sap;
			soil.sompool[SURFFWD].nmass += nmass_litter_sap;
			soil.sompool[SURFFWD].pmass += pmass_litter_sap;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter(cmass_litter_sap, nmass_litter_sap, pmass_litter_sap, SURFFWD);
			}

			if (negligible(soil.sompool[SURFFWD].cmass)) {
				soil.sompool[SURFFWD].ligcfrac = 0.0;
			}
			else {
				double ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SURFFWD].cmass;
				soil.sompool[SURFFWD].ligcfrac = ligcfrac;
			}

			// Fire
			if (firemodel == GLOBFIRM) {
				litterme[SURFFWD] += cmass_litter_sap * pft.pft.litterme;
				fireresist[SURFFWD] += cmass_litter_sap * pft.pft.fireresist;
			}
		}

		// HEART WOOD

		if (!negligible(cmass_litter_heart) || !negligible(nmass_litter_heart) || !negligible(pmass_litter_heart)) {
		//if (!negligible(cmass_litter_heart) || !negligible(nmass_litter_heart)) {

			// Coarse woody debris

			ligcmass_new = cmass_litter_heart * LIGCFRAC_WOOD;
			ligcmass_old = soil.sompool[SURFCWD].cmass * soil.sompool[SURFCWD].ligcfrac;

			// Add to coarse woody pool and update lignin fraction in pool
			soil.sompool[SURFCWD].cmass += cmass_litter_heart;
			soil.sompool[SURFCWD].nmass += nmass_litter_heart;
			soil.sompool[SURFCWD].pmass += pmass_litter_heart;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter(cmass_litter_heart, nmass_litter_heart, pmass_litter_heart, SURFCWD);
			}

			if (negligible(soil.sompool[SURFCWD].cmass)) {
				soil.sompool[SURFCWD].ligcfrac = 0.0;
			}
			else {
				double ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SURFCWD].cmass;
				soil.sompool[SURFCWD].ligcfrac = ligcfrac;
			}

			// Fire
			if (firemodel == GLOBFIRM) {
				litterme[SURFCWD] += cmass_litter_heart * pft.pft.litterme;
				fireresist[SURFCWD] += cmass_litter_heart * pft.pft.fireresist;
			}

		}

		patch.pft.nextobj();
	}

	// FIRE
	if ( firemodel == GLOBFIRM ) {
		if (soil.sompool[SURFSTRUCT].cmass > 0.0) {
			soil.sompool[SURFSTRUCT].litterme   = litterme[SURFSTRUCT]   / soil.sompool[SURFSTRUCT].cmass;
			soil.sompool[SURFSTRUCT].fireresist = fireresist[SURFSTRUCT] / soil.sompool[SURFSTRUCT].cmass;
		}
		if (soil.sompool[SURFMETA].cmass > 0.0) {
			soil.sompool[SURFMETA].litterme   = litterme[SURFMETA]   / soil.sompool[SURFMETA].cmass;
			soil.sompool[SURFMETA].fireresist = fireresist[SURFMETA] / soil.sompool[SURFMETA].cmass;
		}
		if (soil.sompool[SURFFWD].cmass > 0.0) {
			soil.sompool[SURFFWD].litterme    = litterme[SURFFWD]   / soil.sompool[SURFFWD].cmass;
			soil.sompool[SURFFWD].fireresist  = fireresist[SURFFWD] / soil.sompool[SURFFWD].cmass;
		}
		if (soil.sompool[SURFCWD].cmass > 0.0) {
			soil.sompool[SURFCWD].litterme    = litterme[SURFCWD]   / soil.sompool[SURFCWD].cmass;
			soil.sompool[SURFCWD].fireresist  = fireresist[SURFCWD] / soil.sompool[SURFCWD].cmass;
		}
	}

	// Calculate total litter carbon, nitrogen and phosphorus mass for set N:C and P:C ratio of surface microbial pool
	double litter_cmass = soil.sompool[SURFSTRUCT].cmass + soil.sompool[SURFMETA].cmass +
						  soil.sompool[SURFFWD].cmass + soil.sompool[SURFCWD].cmass;
	double litter_nmass = soil.sompool[SURFSTRUCT].nmass + soil.sompool[SURFMETA].nmass +
						  soil.sompool[SURFFWD].nmass + soil.sompool[SURFCWD].nmass;
	double litter_pmass = soil.sompool[SURFSTRUCT].pmass + soil.sompool[SURFMETA].pmass +
						  soil.sompool[SURFFWD].pmass + soil.sompool[SURFCWD].pmass;

	// Set N:C ratio of surface microbial pool based on N:C ratio of litter from all PFTs
	// Parton et al 1993 Fig 4. Dry mass litter == cmass litter * 2
	if (!negligible(litter_cmass)) {
		setntoc(soil, litter_nmass / (litter_cmass * 2.0), SURFMICRO, 20.0, 10.0, 0.0, NCONC_SAT);
		setptoc(soil, litter_pmass / (litter_cmass * 2.0), SURFMICRO, 80.0, 30.0, 0.0, PCONC_SAT);
	}

	// Add litter to solvesom array every month
	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
		if (date.dayofmonth == 0) {
			soil.solvesom.push_back(soil.litterSolveSOM);
			soil.litterSolveSOM.clear();
		}
	}
}


/// LEACHING
/** Leaching fractions for both organic and mineral leaching
 *  Should be called every day in both daily and monthly mode
 */
void leaching(Soil& soil) {

	double minleachfrac;

	if (!negligible(soil.dperc)) {

		// Leaching from available nitrogen mineral pool
		// in proportion to amount of water drainage
		// Use Gerten equivalents here 
		minleachfrac = soil.dperc / (soil.dperc + soil.soiltype.gawc[0] * soil.get_soil_water_upper());

		// Leaching from decayed organic carbon/nitrogen
		// using Parton et al. eqn. 8; CENTURY 5 parameter update; from equation: C Leached=microbial_C*[OMLECH(1)+OMLECH(2)*sand_fraction]*[1.0f-(OMLECH(3)-water_leaching)/OMLECH(3)],
		// reorganised as: leachfrac = [water_leaching/OMLECH(3)]*[OMLECH(1)+OMLECH(2)*sand_fraction]
		// water_leaching was originally in cm/month
		double OMLECH_1 = 0.03;
		double OMLECH_2 = 0.12;
		double OMLECH_3 = 1.9;	// saturation point in leaching equation (cm H2O/month)
		double cmpermonth_to_mmperday = 10.0 * 12.0 / 365.0;

		soil.orgleachfrac = min(1.0, soil.dperc / (OMLECH_3 * cmpermonth_to_mmperday)) * (OMLECH_1 + OMLECH_2 * soil.soiltype.sand_frac);
	}
	else {
		minleachfrac = 0.0;
		soil.orgleachfrac = 0.0;
	}

	// Leaching of soil mineral nitrogen
	// Allowed on days with residual nitrogen following vegetation uptake
	// in proportion to amount of water drainage
	const double nmin_avail = soil.nmass_avail(NO3);
	if (nmin_avail > 0.0) {

		double leaching = nmin_avail * minleachfrac;

		soil.nmass_subtract(leaching,NO3);
		soil.aminleach += leaching;
	}

	// Leaching of soil mineral phosphorus
	// Allowed on days with residual phosphorus following vegetation uptake
	// in proportion to amount of water drainage
	const double pmin_avail = soil.pmass_labile;
	if (pmin_avail > 0.0) {

		double leaching_p = pmin_avail * minleachfrac;

		pmass_add(soil, -leaching_p);

		soil.aminpleach += leaching_p;

	}

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
		soil.morgleach_mean[date.month] += soil.orgleachfrac / date.ndaymonth[date.month];
		soil.mminleach_mean[date.month] += minleachfrac      / date.ndaymonth[date.month];
		soil.mminpleach_mean[date.month] += minleachfrac / date.ndaymonth[date.month];
	}
}


/// Nitrogen addition to the soil
/** Daily nitrogen addition to the soil
 *  from deposition and fixation
 */
void soilnadd(Patch& patch) {

	Soil& soil = patch.soil;

	if (!ifntransform) {
		double nflux = 0.01 * (soil.NH4_input + soil.NO3_input);
		soil.NH4_input -= 0.01 * soil.NH4_input;
		soil.NO3_input -= 0.01 * soil.NO3_input;

		patch.fluxes.report_flux(Fluxes::NH3_SOIL, nflux);

	// Nitrogen deposition and fertilization input to the soil (calculated in snow_ninput())
	}

	// Nitrogen fixation
	// If soil available nitrogen is above the value for minimum SOM C:N ratio, then
	// nitrogen fixation is reduced (nitrogen rich soils)
	const double nmin_avail = soil.nmass_avail(NH4);

	if (nmin_avail < NMASS_SAT) {
		const double daily_nfix = soil.anfix_calc / (double)date.year_length();

		if (nmin_avail + daily_nfix < NMASS_SAT) {
			soil.NH4_mass += daily_nfix;
			soil.anfix += daily_nfix;
		}
		else {
			soil.anfix += NMASS_SAT - soil.NH4_mass;
			soil.NH4_mass += NMASS_SAT - nmin_avail;
		}
	}
	// Nitrogen deposition and fertilization input to the soil (calculated in snow_ninput())
	soil.NH4_mass += soil.NH4_input;
	soil.NO3_mass += soil.NO3_input;

	// Calculate nitrogen fixation (Cleveland et al. 1999)
	// by using five year average aaet
	if (date.islastmonth && date.islastday) {

		// Add this year's AET to aaet_5 which keeps track of the last 5 years
		patch.aaet_5.add(patch.aaet+patch.aevap+patch.aintercep);

		// Calculate estimated nitrogen fixation (aaet should be in cm/yr, eqn is in nitrogen/ha/yr)
		soil.anfix_calc = max((nfix_a * patch.aaet_5.mean() * CM_PER_MM + nfix_b) * HA_PER_M2, 0.0);

		if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
			soil.anfix_mean += soil.anfix;
		}
	}
}

/// Phosphorus addition to the soil
/** Daily phosphorus addition to the soil
*  from deposition and weathering
*/
void soilpadd(Patch& patch) {

	Soil& soil = patch.soil;

	double wfps = soil.wfps(0)*100.0;

	double daily_pwtr, p_temp_effect;

	double R_gas_constant = 8.3144621;
	double ea = patch.get_climate().pwtr_ea;
	double soiltemp = soil.get_soil_temp_25();
	double bi = patch.get_climate().pwtr_bi;
	double pcont = patch.get_climate().pwtr_pcont;
	double shield = patch.get_climate().pwtr_shield;

	if (param["file_pwtr"].str != "") {
		p_temp_effect = exp(-ea / R_gas_constant * (1 / (soiltemp + 273.0) - 1 / 284.15));
		daily_pwtr = bi * (pcont / 100.0) * patch.soil.runoff * p_temp_effect * shield / 1000.0;
	}
	else { 
		daily_pwtr = soil.soiltype.pwtr / date.year_length();
	}


	// Phosphorus weathering input
	pmass_add(soil, daily_pwtr);
	soil.apwtr += daily_pwtr;

	// Phosphorus fertilization and deposition input (calculated in snow_pinput())
	pmass_add(soil, soil.pmass_labile_input);
	soil.apdep += soil.pmass_labile_input;

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
		soil.apwtr_mean += soil.apwtr;
		soil.apdep_mean += soil.apdep;
	}
	


}


/// Vegetation nitrogen uptake
/** Daily vegetation uptake of mineral nitrogen
 *  Partitioned among individuals according to todays nitrogen demand
 *  and individuals root area
 */
void vegetation_n_uptake(Patch& patch) {

	// Daily nitrogen demand given by:
	//	 For individual:
	//     (1)  ndemand = leafndemand + rootndemand + sapndemand;
    //          where
	//          leafndemand is leaf demand based on vmax
	//			rootndemand is based on optimal leaf C:N ratio
	//          sapndemand is based on optimal leaf C:N ratio
	//
	// Actual nitrogen uptake for each day and individual given by:
	//     (2)  nuptake = ndemand * fnuptake
	//     where fnuptake is individual uptake capacity calculated in fnuptake in canexch.cpp

	double nuptake_day;
	double nuptake_day_myco;
	double nmass_herb;

	Vegetation& vegetation=patch.vegetation;
	Soil& soil = patch.soil;

	//const double orignmass = soil.NH4_mass + soil.NO3_mass;
	const double orignmass = soil.nmass_avail();
	double ammonium_frac = orignmass ? soil.NH4_mass / orignmass : 0;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (date.day == 0)
			indiv.anuptake = 0.0;

		nuptake_day_myco = 0.0;

		if (indiv.myco_type)
			//nuptake_day_myco = indiv.ndemand * indiv.fractomax_nmyco;
			nuptake_day_myco = indiv.ndemand * indiv.fnuptake * indiv.fractomax_nmyco;

		nuptake_day           = indiv.ndemand * indiv.fnuptake;
		//nuptake_day = indiv.ndemand * indiv.fnuptake + nuptake_day_myco;

		if (ifherbivory)
			nmass_herb = indiv.leaffndemand  * nuptake_day * indiv.herb_frac;
		else 
			nmass_herb = 0.0;

		indiv.anuptake        += nuptake_day;
		indiv.nmass_leaf      += indiv.leaffndemand  * nuptake_day - nmass_herb;
		indiv.nmass_root      += indiv.rootfndemand  * nuptake_day;
		indiv.nmass_sap       += indiv.sapfndemand   * nuptake_day;
		if (indiv.pft.phenology == CROPGREEN && ifnlim)
			indiv.cropindiv->nmass_agpool += indiv.storefndemand * nuptake_day;
		else
			indiv.nstore_longterm += indiv.storefndemand * nuptake_day;


		if (ifntransform) {
			if (!indiv.myco_type) {
				double ammonium = (nuptake_day - nmass_herb) * ammonium_frac;  // Only no errors if ammonium frac is considered...
				soil.NH4_mass -= ammonium;
				soil.NO3_mass -= (nuptake_day - nmass_herb) - ammonium;
				//double ammonium = (nuptake_day - nmass_herb) * indiv.frac_nh4;  // Only no errors if ammonium frac is considered...
				//soil.NH4_mass -= ammonium;
				//soil.NO3_mass -= (nuptake_day - nmass_herb) - ammonium;
			}
			else {
				double ammonium = (nuptake_day - nuptake_day_myco - nmass_herb) * ammonium_frac;  // Only no errors if ammonium frac is considered...
				soil.NH4_mass -= ammonium;
				soil.NO3_mass -= (nuptake_day - nuptake_day_myco - nmass_herb) - ammonium;
				//double ammonium = (nuptake_day - nuptake_day_myco - nmass_herb) * indiv.frac_nh4;  // Only no errors if ammonium frac is considered...
				//soil.NH4_mass -= ammonium;
				//soil.NO3_mass -= (nuptake_day - nuptake_day_myco - nmass_herb) - ammonium;
				/*soil.sompool[SOILSTRUCT].nmass -= nuptake_day_myco * 0.18;
				soil.sompool[SOILMETA].nmass -= nuptake_day_myco * 0.82;*/
				soil.sompool[SURFSTRUCT].nmass -= nuptake_day_myco * 0.18;
				soil.sompool[SURFMETA].nmass -= nuptake_day_myco * 0.82;
				/*soil.sompool[SURFSTRUCT].nmass -= nuptake_day_myco * 0.10;
				soil.sompool[SURFMETA].nmass -= nuptake_day_myco * 0.20;
				soil.sompool[SURFCWD].nmass -= nuptake_day_myco * 0.35;
				soil.sompool[SURFFWD].nmass -= nuptake_day_myco * 0.35;*/

			}

		} 
		else {
			soil.nmass_subtract(nuptake_day - nuptake_day_myco - nmass_herb);
		}

		// Report herbivory nitrogen flux
		indiv.report_flux(Fluxes::N_HERB, nmass_herb);

		/*soil.sompool[SOILSTRUCT].nmass -= nuptake_day_myco * 0.18;
		soil.sompool[SOILMETA].nmass -= nuptake_day_myco * 0.82;*/

		if (!negligible(indiv.phen))
			indiv.cton_leaf_aavr += min(indiv.cton_leaf(),indiv.cton_leaf_max);

		vegetation.nextobj();
	}

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr && !negligible(orignmass)) {
		soil.fnuptake_mean[date.month] += (1.0 - soil.nmass_avail() / orignmass) / date.ndaymonth[date.month];
	}
}

/// Vegetation phosphorus uptake
/** Daily vegetation uptake of mineral phosphorus
*  Partitioned among individuals according to todays phosphorus demand
*  and individuals root area
*/
void vegetation_p_uptake(Patch& patch) {

	// Daily phosphorus demand given by:
	//	 For individual:
	//     (1)  pdemand = leafpdemand + rootpdemand + sappdemand;
	//          where
	//          leafpdemand is leaf demand based on vmax
	//			rootpdemand is based on optimal leaf C:P ratio
	//          sapnpemand is based on optimal leaf C:P ratio
	//
	// Actual phosphorus uptake for each day and individual given by:
	//     (2)  puptake = pdemand * fpuptake
	//     where fpuptake is individual uptake capacity calculated in fpuptake in canexch.cpp

	double puptake_day;
	double puptake_day_myco;
	double pmass_herb;

	Vegetation& vegetation = patch.vegetation;
	Soil& soil = patch.soil;

	const double origpmass = soil.pmass_labile;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (date.day == 0)
			indiv.apuptake = 0.0;

		puptake_day_myco = 0.0;

		if (indiv.myco_type)
			//puptake_day_myco = indiv.pdemand * indiv.fractomax_pmyco;
			puptake_day_myco = indiv.pdemand * indiv.fpuptake * indiv.fractomax_pmyco;

		//puptake_day = indiv.pdemand * indiv.fpuptake + puptake_day_myco;
		puptake_day = indiv.pdemand * indiv.fpuptake;

		if (ifherbivory)
			pmass_herb = indiv.leaffndemand  * puptake_day * indiv.herb_frac;
		else 
			pmass_herb = 0.0;

		indiv.apuptake += puptake_day;
		indiv.pmass_leaf += indiv.leaffpdemand  * puptake_day - pmass_herb;
		indiv.pmass_root += indiv.rootfpdemand  * puptake_day;
		indiv.pmass_sap += indiv.sapfpdemand   * puptake_day;
		if (indiv.pft.phenology == CROPGREEN && ifplim)
			indiv.cropindiv->pmass_agpool += indiv.storefpdemand * puptake_day;
		else
			indiv.pstore_longterm += indiv.storefpdemand * puptake_day;
		
		//soil.pmass_labile_delta -= (puptake_day - puptake_day_myco - pmass_herb);
		pmass_add(soil, -puptake_day + puptake_day_myco + pmass_herb);

		/*soil.sompool[SOILSTRUCT].pmass -= puptake_day_myco * 0.18;
		soil.sompool[SOILMETA].pmass -= puptake_day_myco * 0.82;*/
		soil.sompool[SURFSTRUCT].pmass -= puptake_day_myco * 0.18;
		soil.sompool[SURFMETA].pmass -= puptake_day_myco * 0.82;
		/*soil.sompool[SURFSTRUCT].pmass -= puptake_day_myco * 0.10;
		soil.sompool[SURFMETA].pmass -= puptake_day_myco * 0.20;
		soil.sompool[SURFCWD].pmass -= puptake_day_myco * 0.35;
		soil.sompool[SURFFWD].pmass -= puptake_day_myco * 0.35;*/
				
		if (!negligible(indiv.phen))
			indiv.ctop_leaf_aavr += min(indiv.ctop_leaf(), indiv.ctop_leaf_max);

		// Report herbivory phosphorus flux
		indiv.report_flux(Fluxes::P_HERB, pmass_herb);

		vegetation.nextobj();
	}

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr && !negligible(origpmass)) {
		soil.fpuptake_mean[date.month] += (1.0 - soil.pmass_labile / origpmass) / date.ndaymonth[date.month];
	}
}


/// Add litter to equilsom()
void add_litter(Soil& soil, int year, int pool) {

	// If LUC have occurred in between solvesomcent_beginyr and solvesomcent_endyr then array might be too smalll
	if (year >= (int)soil.solvesom.size())
		year = year%(int)soil.solvesom.size();

	// Add to litter pools
	soil.sompool[pool].cmass += soil.solvesom[year].get_clitter(pool);
	soil.sompool[pool].nmass += soil.solvesom[year].get_nlitter(pool);
	soil.sompool[pool].pmass += soil.solvesom[year].get_plitter(pool);
}

/// Iteratively solving differential flux equations for century SOM pools
/** Iteratively solving differential flux equations for century SOM pools
 *  assuming annual litter inputs, nitrogen uptake and leaching is close
 *  to long term equilibrium
 */
void equilsom(Soil& soil) {

	if (!(date.year == soil.solvesomcent_endyr && date.islastmonth && date.islastday))
		return;

	// Number of years to run SOM pools, value chosen to get cold climates to equilibrium
	const int EQUILSOM_YEARS = 40000;

	Patch& patch = soil.patch;
	const Climate& climate = soil.patch.get_climate();
	const Gridcell& gridcell = patch.stand.get_gridcell();

	// Save soil mineral nitrogen status
	double save_NH4_mass = soil.nmass_avail(NH4);
	double save_NO3_mass = 0.0;
	if (ifntransform) {
		save_NO3_mass = soil.nmass_avail(NO3);
	}

	// Save pmass_labile status
	double save_pmass_labile = soil.pmass_labile;

	// Number of years with mean input data
	int nyear = soil.solvesomcent_endyr - soil.solvesomcent_beginyr + 1;

	for (int m = 0; m < 12; m++) {

		// Monthly average decay rates
		for (int p = 0; p < NSOMPOOL; p++) {
			soil.sompool[p].mfracremain_mean[m] = pow(soil.sompool[p].mfracremain_mean[m] / nyear, (double)date.ndaymonth[m]);
		}

		// Monthly average mineral nitrogen uptake
		soil.fnuptake_mean[m] /= nyear;

		// Monthly average mineral phosphorus uptake
		soil.fpuptake_mean[m] /= nyear;

		// Monthly average organic carbon and nitrogen leaching
		soil.morgleach_mean[m] /= nyear;

		// Monthly average mineral nitrogen leaching
		soil.mminleach_mean[m] /= nyear;

		// Monthly average mineral phosphorus leaching
		soil.mminpleach_mean[m] /= nyear;
	}

	// Annual average nitrogen fixation
	soil.anfix_mean /= nyear;

	// Annual phosphorus weathering
	soil.apwtr_mean /= nyear;

	// Spin SOM pools with saved litter input, nitrogen addition and fractions of
	// nitrogen uptake and leaching for EQUILSOM_YEARS years with monthly timesteps
	for (int yr = 0; yr < EQUILSOM_YEARS; yr++) {

		// Which year in saved data set
		int savedyear = yr%nyear;

		// Monthly time steps
		for (int m = 0; m < 12; m++) {

			// Transfer yearly mean litter on first day of year
			add_litter(soil, savedyear*12+m, SURFSTRUCT);
			add_litter(soil, savedyear*12+m, SURFMETA);
			add_litter(soil, savedyear*12+m, SOILSTRUCT);
			add_litter(soil, savedyear*12+m, SOILMETA);
			add_litter(soil, savedyear*12+m, SURFFWD);
			add_litter(soil, savedyear*12+m, SURFCWD);

			// Calculate total litter carbon and nitrogen mass for set N:C ratio of surface microbial pool
			double litter_cmass = soil.sompool[SURFSTRUCT].cmass + soil.sompool[SURFMETA].cmass +
				soil.sompool[SURFFWD].cmass + soil.sompool[SURFCWD].cmass;
			double litter_nmass = soil.sompool[SURFSTRUCT].nmass + soil.sompool[SURFMETA].nmass +
				soil.sompool[SURFFWD].nmass + soil.sompool[SURFCWD].nmass;
			double litter_pmass = soil.sompool[SURFSTRUCT].pmass + soil.sompool[SURFMETA].pmass +
				soil.sompool[SURFFWD].pmass + soil.sompool[SURFCWD].pmass;

			// Set N:C ratio of surface microbial pool based on N:C ratio of litter from all PFTs
			if (!negligible(litter_cmass)) {
				setntoc(soil, litter_nmass / (litter_cmass * 2.0), SURFMICRO, 20.0, 10.0, 0.0, NCONC_SAT);
				setptoc(soil, litter_pmass / (litter_cmass * 2.0), SURFMICRO, 80.0, 30.0, 0.0, PCONC_SAT);
			}

			// Monthly nitrogen uptake
			soil.NH4_mass *= (1.0 - soil.fnuptake_mean[m]);
			soil.NO3_mass *= (1.0 - soil.fnuptake_mean[m]);

			// Monthly mineral nitrogen leaching
			double nmass = soil.nmass_avail(NO3);
			double nleach = nmass * (1.0 - soil.mminleach_mean[m]);
			soil.nmass_subtract(nleach, NO3);

			// Monthly phosphorus uptake
			soil.pmass_labile *= (1.0 - soil.fpuptake_mean[m]);

			// Monthly mineral phosphorus leaching
			 double pleach = soil.pmass_labile * (1.0 - soil.mminpleach_mean[m]);
			soil.pmass_labile -= pleach;
			if (soil.pmass_labile < 0.0) {
				patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_labile);
				soil.pmass_labile = 0.0;
			}

			// Monthly nitrogen addition to the system
			soil.nmass_inc((gridcell.aNH4dep + soil.anfix_mean) / 12.0, NH4);
			soil.nmass_inc(gridcell.aNO3dep / 12.0, NO3);

			// Monthly phosphorus addition to the system
			soil.pmass_labile += (soil.apwtr_mean + soil.apdep_mean) / 12.0;

			// Monthly decomposition and fluxes between SOM pools

			// Set this months decay rates
			for (int p = 0; p < NSOMPOOL; p++) {
				soil.sompool[p].fracremain = soil.sompool[p].mfracremain_mean[m];
			}

			// Set this months organic nitrogen leaching fraction
			soil.orgleachfrac = soil.morgleach_mean[m];

			// Monthly decomposition and fluxes between SOM pools
			// and nitrogen flux from soil
			somfluxes(patch, true, false);
		}
	}

	// Reset mineral nitrogen status
	soil.NH4_mass = save_NH4_mass;
	soil.NO3_mass = save_NO3_mass;

	// Reset pmass_labile status
	soil.pmass_labile = save_pmass_labile;

	// Reset variables for next equilsom()
	for (int m = 0; m < 12; m++) {

		for (int p = 0; p < NSOMPOOL; p++) {
			soil.sompool[p].mfracremain_mean[m] = 0.0;
		}

		soil.fnuptake_mean[m]  = 0.0;
		soil.morgleach_mean[m] = 0.0;
		soil.mminleach_mean[m] = 0.0;
		soil.mminpleach_mean[m] = 0.0;
	}

	soil.anfix_mean = 0.0;
	soil.apwtr_mean = 0.0;
	soil.apdep_mean = 0.0;
	soil.solvesom.clear();
	std::vector<LitterSolveSOM>().swap(soil.solvesom); // clear array memory
}

/// SOM CENTURY DYNAMICS
/** To be called each simulation day for each modelled patch, following update
 *  of soil temperature and soil water.
 *  Transfers litter, performes nitrogen uptake and addition, leaching and decomposition.
 */
void som_dynamics_century(Patch& patch, Climate& climate, bool tillage) {

	// Transfer litter to SOM pools
	transfer_litter(patch);

	// Daily nitrogen uptake
	vegetation_n_uptake(patch);

	// Daily phosphorus uptake
	vegetation_p_uptake(patch);

	// Daily nitrogen addition to the soil
	soilnadd(patch);

	// Daily phosphorus addition to the soil
	soilpadd(patch);

	// Daily mineral and organic nitrogen leaching
	leaching(patch.soil);

	// Daily or monthly decomposition and fluxes between SOM pools
	somfluxes(patch, false, tillage);

	// Nitrogen transformation in soil
	ntransform(patch, climate);

	// Solve SOM pool sizes at end of year given by soil.solvesomcent_endyr
	equilsom(patch.soil);
}

/// Choose between CENTURY or standard LPJ SOM dynamics
/**
*/
void som_dynamics(Patch& patch, Climate& climate) {

	bool tillage = iftillage && patch.stand.landcover == CROPLAND;
	if (ifcentury) {
		som_dynamics_century(patch, climate, tillage);
	}
	else {
		som_dynamics_lpj(patch, tillage);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Chatskikh, D., Hansen, S., Olesen, J.E. & Petersen, B.M. 2009. A simplified modelling approach
//	 for quantifying tillage effects on soil carbon stocks. Eur.J.Soil.Sci., 60:924-934.
// CENTURY Soil Organic Matter Model Version 5; Century User's Guide and Reference.
//  http://www.nrel.colostate.edu/projects/century5/reference/index.htm; accessed Oct.7, 2016.
// Comins, H. N. & McMurtrie, R. E. 1993. Long-Term Response of Nutrient-Limited
//   Forests to CO2 Enrichment - Equilibrium Behavior of Plant-Soil Models.
//   Ecological Applications, 3, 666-681.
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical exploration
//   of the relationships of soil moisture characteristic to the physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Cleveland C C (1999) Global patterns of terrestrial biological nitrogen (N2) fixation
//   in natural ecosystems. GBC 13: 623-645
// Foley J A 1995 An equilibrium model of the terrestrial carbon budget
//   Tellus (1995), 47B, 310-319
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A
//   process-based, terrestrial biosphere model of ecosystem dynamics
//   (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Kirschbaum, M. U. F. and K. I. Paul (2002). "Modelling C and N dynamics in forest soils
//   with a modified version of the CENTURY model." Soil Biology & Biochemistry 34(3): 341-354.
// Meentemeyer, V. (1978) Macroclimate and lignin control of litter decomposition
//   rates. Ecology 59: 465-472.
// Mikan, C. J., Schimel, J. P., and Doyle, A. P. 2002. Temperature controls of microbial respiration 
//   in arctic tundra soils above and below freezing, Soil Biol.Biochem., 34, 1785–1795, 2002.
// Parton, W. J., Scurlock, J. M. O., Ojima, D. S., Gilmanov, T. G., Scholes, R. J., Schimel, D. S.,
//   Kirchner, T., Menaut, J. C., Seastedt, T., Moya, E. G., Kamnalrut, A. & Kinyamario, J. I. 1993.
//   Observations and Modeling of Biomass and Soil Organic-Matter Dynamics for the Grassland Biome
//   Worldwide. Global Biogeochemical Cycles, 7, 785-809.
// Parton, W. J., Hanson, P. J., Swanston, C., Torn, M., Trumbore, S. E., Riley, W. & Kelly, R. 2010.
//   ForCent model development and testing using the Enriched Background Isotope Study experiment.
//   Journal of Geophysical Research-Biogeosciences, 115.
// Schaefer, K. and Jafarov, E. 2016. 
//   A parameterization of respiration in frozen soils based on substrate availability
//   Biogeosciences, 13, 1991 - 2001, https ://doi.org/10.5194/bg-13-1991-2016

