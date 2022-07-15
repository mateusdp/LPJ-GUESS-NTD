////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file management.cpp
/// \brief Management functions for cropland, managed forest and pasture
/// \author Mats Lindeskog
/// $Date:  $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "landcover.h"
#include "management.h"
#include "driver.h"

// Automated clearcut density targets for needle-leaf and broad-leaf trees. These values are sensitive to changes in 
// PFT growth and mortality and should be calibrated for each PFT and geographic region for realistic rotation times.
const int DENSTARGET_NL = 400;	// Bellassen (2010) value 100, modified for LPJ-GUESS simulations.
const int DENSTARGET_BL = 100;	// Bellassen (2010) value 200, modified for LPJ-GUESS simulations.

// Functions to check available wood for harvest at individual, patch and stand levels

/// Checks available wood for harvest at individual level.
/** With default parameters, the function returns total potential harvested C for an individual, both potential
 *  harvest products and potential fuel wood+residues. Harvest efficiency is accounted for.
 *  
 *  INPUT PARAMETERS
 *  \param stem_cmass_only				if true, function returns potential removed stem C
 *  \param to_product_pool				if true, function returns potential removed stem products C
 */
double check_harvest_cmass(Individual& indiv, bool stem_cmass_only, bool to_product_pool) {

	double cmass_harvest;
	Stand& stand = indiv.vegetation.patch.stand;
	Harvest_CN cp;

	cp.copy_from_indiv(indiv, indiv.has_daily_turnover(), false);

	// Default forestry harvest parameters: pft values
	double harv_eff_wood_harvest = indiv.pft.harv_eff;				// 0.9
	double res_outtake_twig_wood_harvest = indiv.pft.res_outtake;	// 0.4
	double res_outtake_coarse_root_wood_harvest = 0.1;
	// Use ManagementType values, if defined
	ManagementType& mt = indiv.vegetation.patch.stand.get_current_management();
	if(mt.harv_eff_cc != -1.0)
		harv_eff_wood_harvest = mt.harv_eff_cc;
	if(mt.res_outtake_twig_cc != -1.0)
		res_outtake_twig_wood_harvest = mt.res_outtake_twig_cc;
	if(mt.res_outtake_coarse_root_cc != -1.0)
		res_outtake_coarse_root_wood_harvest = mt.res_outtake_coarse_root_cc;

	// Harvest of transferred areas:
	harvest_wood(cp, indiv.diam, indiv.pft, indiv.alive, 1.0, harv_eff_wood_harvest, res_outtake_twig_wood_harvest,
					res_outtake_coarse_root_wood_harvest);
	if(to_product_pool) {
		cmass_harvest = cp.acflux_harvest_wood_toprod;
	}
	else if(stem_cmass_only) {
		cmass_harvest = cp.acflux_harvest_wood;
	}
	else {
		cmass_harvest = (cp.acflux_harvest + cp.acflux_harvest_wood_toprod);
	}

	return cmass_harvest;
}

/// Checks available wood for harvest at patch level.
/** With default parameters, the function returns total potential harvested C for a patch, both potential harvest
 *  products and potential fuel wood+residues. Harvest efficiency is accounted for.
 *  
 *  INPUT PARAMETERS
 *  \param stem_cmass_only				if true, function returns potential removed stem C
 *  \param check_selection				if true, function check only value for PFTs in selection
 */
double check_harvest_cmass(Patch& patch, bool stem_cmass_only, bool check_selection) {

	double cmass_harvest = 0.0;
	ManagementType& mt = patch.stand.get_current_management();
	StandType& st = stlist[patch.stand.stid];
	Vegetation& vegetation = patch.vegetation;

	for(unsigned int i=0;i<vegetation.nobj;i++) {
		Individual& indiv = vegetation[i];
		Standpft& spft = patch.stand.pft[indiv.pft.id];

		if(!check_selection || spft.plant || mt.planting_system == "")
			cmass_harvest += check_harvest_cmass(indiv, stem_cmass_only);
	}
	return cmass_harvest;
}

/// Checks available wood for harvest at stand level.
/** With default parameters, the function returns total potential harvested C for a stand, both potential harvest
 *  products and potential fuel wood+residues. Harvest efficiency is accounted for.
 *  
 *  INPUT PARAMETERS
 *  \param stem_cmass_only				if true, function returns potential removed stem C
 *  \param check_selection				if true, check only value for PFTs in selection
 */
double check_harvest_cmass(Stand& stand, bool stem_cmass_only, bool check_selection) {

	double cmass_harvest = 0.0;

	for(unsigned int i=0;i<stand.nobj;i++) {
		Patch& patch = stand[i];
		cmass_harvest += check_harvest_cmass(patch, stem_cmass_only, check_selection) / (double)stand.nobj;
	}
	return cmass_harvest;
}

/// Harvest function used for managed forest and for clearing natural vegetation during land use change.
/** This version does not directly change the biomass and litter in an individual and the corresponding patch/patchpft,
 *  but updates values for a Harvest_CN struct copy.
 *  A fraction of trees is cut down (frac_cut)
 *  A fraction of stem wood is harvested (pft.harv_eff). A fraction of harvested wood (pft.harvest_slow_frac or 100% of
 *  trees above a diameter limit if harvest_burn_thin_trees == true) is sent to harvested_products_slow and the rest plus
 *  residue outtake is sent to acflux_harvest.The rest, including leaves and roots, is sent to litter.
 *  Called from landcover_dynamics() first day of the year if any natural or forest vegetation is transferred to another,
 *  land use or from harvest_wood(Individual& indiv, ...).
 *
 *  INPUT PARAMETERS
 *  \param diam						tree diameter (m)
 *  \param frac_cut					fraction of trees cut
 *  \param harv_eff					harvest efficiency
 *  \param res_outtake_twig			removed twig fraction
 *  \param res_outtake_coarse_root	removed course root fraction
 *									global parameters:
 *   - harvest_burn_thin_trees		whether to put all harvested wood from trees with diam > a diameter limit in product
 *									pool and burn all wood of trees below the same diameter limit.
 *   - ifslowharvestpool			whether a fraction of harvested wood is put into a product pool
 *  \param pft						reference to a Pft containing the following public members:
 *   - stem_frac    				fraction of wood cmass that belongs to stems
 *   - twig_frac    				fraction of wood cmass that belongs to twigs
 *   - lifeform	    				life form (tree or grass)
 *   - harvest_slow_frac			fraction of harvested products that goes to long-lived products
 *   - leafphysiognomy 				leaf physiognomy (needleleaf, broadleaf)
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members copied to and from the corresponding
 *									variables of an Individual:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgN/m2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgN/m2)
 *   - nstore_labile    			labile nitrogen storage (kgN/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgN/m2)
 *   - max_n_storage				maximum size of nitrogen storage
 *									,and the following public members copied to and from the corresponding variables of
 *									a Patchpft containing an Individual:
 *   - litter_leaf    				leaf C litter (kgC/m2)
 *   - litter_root 					root C litter (kgC/m2)
 *   - litter_sap   				sapwood C litter (kgC/m2)
 *   - litter_heart   				heartwood C litter (kgC/m2)
 *   - nmass_litter_leaf 			leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			root nitrogen litter (kgN/m2)
 *   - nmass_litter_sap 			sapwood nitrogen litter (kgN/m2)
 *   - nmass_litter_heart        	heartwood nitrogen litter (kgN/m2)
 *   - cmass_harvested_products_slow wood product pool (kgC/m2)
 *   - nmass_harvested_products_slow wood product pool nitrogen (kgN/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members added to the corresponding variables of
 *									a Patchpft corresponding to an Individual:
 *   - acflux_harvest_wood			harvested wood C before removing part to the product pool (kgC/m2)
 *   - acflux_harvest_wood_toprod	harvested wood C removed to the product pool (kgC/m2)
 *   - acflux_harvest_tolitter		harvested tree C left as litter (kgC/m2)
 *									,and the following public members members added to the corresponding variables of a
 *									Patch containing an individual::
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgN/m2)
 */
void harvest_wood(Harvest_CN& i, double diam, Pft& pft, bool alive, double frac_cut, double harv_eff, double res_outtake_twig,
					double res_outtake_coarse_root) {

	// only harvest trees
	if (pft.lifeform == GRASS)
		return;

	double stem_harvest = 0.0;
	double residue_outtake = 0.0;
	// Fraction of wood cmass that are stems
	double stem_frac = pft.stem_frac; // Default 0.65
	// Fraction of wood cmass that are twigs
	double twig_frac = pft.twig_frac; // Default 0.13
	// Fraction of wood cmass that are coarse roots, including stumps
	double coarse_root_frac = 1.0 - stem_frac - twig_frac;	// 0.22 with default stem_frac and twig_frac values
	// Fraction of leaves adhering to twigs at the time of removal (provisional value, typical for situation for harvested spruce in Sweden)
	double adhering_leaf_frac = 0.75;

	double harvest_slow_frac = pft.harvest_slow_frac;	// Default 0.33

	if(harvest_burn_thin_trees) {
		// Put all harvested wood into product pool above a diameter limit, burn the smaller trees (limits from Anne-Sofie
		// Lanso, personal communication, as used with ORCHIDEE, except for tropical broadleaves, which is the reported
		// limit for general harvesting purposes in Guyana (R.Heinrich, FAO: SUSTAINABLE FOREST HARVESTING 1) 
		// Whether the pft is boreal, temperate or tropical is determined by the pft.pstemp_low parameter value.
		const double DIAMETER_LIMIT_NEEDLELAF = 0.2;
		const double DIAMETER_LIMIT_BROADLEAF_BOREAL = 0.2;
		const double DIAMETER_LIMIT_BROADLEAF_TEMPERATE = 0.3;
		const double DIAMETER_LIMIT_BROADLEAF_TROPICAL = 0.35;
		if(pft.leafphysiognomy == NEEDLELEAF && diam < DIAMETER_LIMIT_NEEDLELAF ||
			pft.leafphysiognomy == BROADLEAF && pft.is_boreal() && diam < DIAMETER_LIMIT_BROADLEAF_BOREAL ||
			pft.leafphysiognomy == BROADLEAF && pft.is_temperate() && diam < DIAMETER_LIMIT_BROADLEAF_TEMPERATE ||
			pft.leafphysiognomy == BROADLEAF && pft.is_tropical() && diam < DIAMETER_LIMIT_BROADLEAF_TROPICAL ) {
			harvest_slow_frac = 0.0;
		}
		else {
			harvest_slow_frac = 1.0;
		}
	}

	// all root carbon and nitrogen goes to litter
	if (alive) {

		i.cmass_litter_root += i.cmass_root * frac_cut;
		i.acflux_harvest_tolitter += i.cmass_root * frac_cut;
		i.cmass_root *= (1.0 - frac_cut);
	}

	i.nmass_litter_root += i.nmass_root * frac_cut;
	i.nmass_litter_root += (i.nstore_labile + i.nstore_longterm) * frac_cut;
	i.nmass_root *= (1.0 - frac_cut);
	i.nstore_labile *= (1.0 - frac_cut);
	i.nstore_longterm *= (1.0 - frac_cut);

	if (alive) {

		// Carbon:

		if (i.cmass_debt <= i.cmass_sap + i.cmass_heart) {

			// harvested stem wood
			stem_harvest += harv_eff * stem_frac * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;
			i.acflux_harvest_wood += stem_harvest;

			// harvested products not consumed (oxidised) this year put into harvested_products_slow
			if (ifslowharvestpool) {
				i.cmass_harvested_products_slow += stem_harvest * harvest_slow_frac;
				i.acflux_harvest_wood_toprod += stem_harvest * harvest_slow_frac;
				stem_harvest *= (1.0 - harvest_slow_frac);
			}

			// harvested products consumed (oxidised) this year put into acflux_harvest
			i.acflux_harvest += stem_harvest;

			// removed leaves adhering to twigs
			residue_outtake += res_outtake_twig * adhering_leaf_frac * i.cmass_leaf * frac_cut;

			// removed twigs
			residue_outtake += res_outtake_twig * twig_frac * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;

			// removed coarse roots
			residue_outtake += res_outtake_coarse_root * coarse_root_frac * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;

			// removed residues are oxidised
			i.acflux_harvest += residue_outtake;

			// not removed residues are put into litter
			i.cmass_litter_leaf += i.cmass_leaf * (1.0 - res_outtake_twig * adhering_leaf_frac) * frac_cut;
			i.acflux_harvest_tolitter += i.cmass_leaf * (1.0 - res_outtake_twig * adhering_leaf_frac) * frac_cut;

			double to_partition_sap   = 0.0;
			double to_partition_heart = 0.0;

			if (i.cmass_heart >= i.cmass_debt) {
				to_partition_sap   = i.cmass_sap;
				to_partition_heart = i.cmass_heart - i.cmass_debt;
			}
			else {
				to_partition_sap   = i.cmass_sap + i.cmass_heart - i.cmass_debt;
			//	dprintf("ATTENTION: pft %s: cmass_debt > cmass_heart; difference=%f\n", (char*)pft.name, i.cmass_debt-i.cmass_heart);
			}
			double wood_residue_frac_to_litter = (1.0 - res_outtake_twig * twig_frac - res_outtake_coarse_root * coarse_root_frac
				- harv_eff * stem_frac) * frac_cut;
			i.cmass_litter_sap += to_partition_sap * wood_residue_frac_to_litter;
			i.acflux_harvest_tolitter += to_partition_sap * wood_residue_frac_to_litter;
			i.cmass_litter_heart += to_partition_heart * wood_residue_frac_to_litter;
			i.acflux_harvest_tolitter += to_partition_heart * wood_residue_frac_to_litter;
		}
		// debt larger than existing wood biomass
		else {
			double debt_excess = i.cmass_debt - (i.cmass_sap + i.cmass_heart);
			dprintf("ATTENTION: cmass_debt > i.cmass_sap + i.cmass_heart; debt_excess=%f\n", debt_excess);
		//	i.debt_excess += debt_excess * frac_cut;	// debt_excess currently not dealt with during wood harvest
		}

		// unharvested trees:
		i.cmass_leaf *= (1.0 - frac_cut);
		i.cmass_sap *= (1.0 - frac_cut);
		i.cmass_heart *= (1.0 - frac_cut);
		i.cmass_debt *= (1.0 - frac_cut);

		//Nitrogen:

		stem_harvest = 0.0;

		// harvested products
		stem_harvest += harv_eff * stem_frac * (i.nmass_sap + i.nmass_heart) * frac_cut;

		// harvested products not consumed this year put into harvested_products_slow_nmass
		if (ifslowharvestpool) {
			i.nmass_harvested_products_slow += stem_harvest * harvest_slow_frac;
			stem_harvest = stem_harvest * (1.0 - harvest_slow_frac);
		}

		// harvested products consumed this year put into anflux_harvest
		i.anflux_harvest += stem_harvest;

		residue_outtake = 0.0;

		// removed leaves adhering to twigs
		residue_outtake += res_outtake_twig * adhering_leaf_frac * i.nmass_leaf * frac_cut;

		// removed twigs
		residue_outtake += res_outtake_twig * twig_frac * (i.nmass_sap + i.nmass_heart) * frac_cut;

		// removed coarse roots
		residue_outtake += res_outtake_coarse_root * coarse_root_frac * (i.nmass_sap + i.nmass_heart) * frac_cut;

		// removed residues are oxidised
		i.anflux_harvest += residue_outtake;

		// not removed residues are put into litter
		i.nmass_litter_leaf += i.nmass_leaf * (1.0 - res_outtake_twig * adhering_leaf_frac) * frac_cut;
		i.nmass_litter_sap += i.nmass_sap * (1.0 - res_outtake_twig * twig_frac - res_outtake_coarse_root * coarse_root_frac
			- harv_eff * stem_frac) * frac_cut;
		i.nmass_litter_heart += i.nmass_heart * (1.0 - res_outtake_twig * twig_frac - res_outtake_coarse_root * coarse_root_frac
			- harv_eff * stem_frac) * frac_cut;

		// unharvested trees:
		i.nmass_leaf *= (1.0 - frac_cut);
		i.nmass_sap *= (1.0 - frac_cut);
		i.nmass_heart *= (1.0 - frac_cut);

		// rescale max_n_storage
		i.max_n_storage *= (1.0 - frac_cut);
	}
}

/// Harvest function used for managed forest and for clearing natural vegetation at land use change
/** This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_wood function.
 *  After the execution of the main harvest_wood function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  A fraction of trees is cut down (frac_cut)
 *  A fraction of stem wood is harvested (pft.harv_eff). A fraction of harvested wood (pft.harvest_slow_frac or 100% of
 *  trees above a diameter limit if harvest_burn_thin_trees == true) is returned as harvested_products_slow and the rest
 *  plus residue outtake is returned as acflux_harvest.The rest, including leaves and roots, is returned as litter.
 *  Called from harvest_forest() first day of the year for normal wood harvest.
 *  Also called first day of the year from set_management() if natural or forest stand is transferred to another management
 *  by cloning (LUC) or after management rotation.
 *  If lc_change is true, harvest C fluxes go to lc.acflux_wood_harvest, if false, they go to Fluxes::HARVESTC.
 *
 *  INPUT PARAMETERS
 *  \param frac_cut					fraction of trees cut
 *  \param harv_eff					harvest efficiency
 *  \param res_outtake_twig			removed twig fraction
 *  \param res_outtake_coarse_root	removed course root fraction
 *  \param lc_change				whether to save harvest in gridcell-level lc struct
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgC/m2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgC/m2)
 *   - nstore_labile    			labile nitrogen storage (kgC/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgC/m2)
 *									Patchpft (accessed from indiv) public members:
 *   - litter_leaf    				leaf C litter (kgC/m2)
 *   - litter_root 					root C litter (kgC/m2)
 *   - litter_sap   				sapwood C litter (kgC/m2)
 *   - litter_heart   				heartwood C litter (kgC/m2)
 *   - nmass_litter_leaf 			leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			root nitrogen litter (kgN/m2)
 *   - nmass_litter_sap 			sapwood nitrogen litter (kgN/m2)
 *   - nmass_litter_heart        	heartwood nitrogen litter (kgN/m2)
 *   - harvested_products_slow		wood product pool (kgC/m2)
 *   - harvested_products_slow_nmass wood product pool nitrogen (kgN/m2)
 *   - cmass_wood_harv        		harvested wood C before removing part to the product pool (kgC/m2)
 *   - cmass_wood_harv_toprod       harvested wood C removed to the product pool (kgC/m2)
 *   - cmass_harv_tolitter        	harvested tree C left as litter (kgC/m2)
 *									Patch (accessed from indiv) public members:
 *   - Fluxes::HARVESTC				harvest flux to atmosphere (kgC/m2)
 *   - Fluxes::HARVESTN   			harvest nitrogen flux out of system (kgN/m2)
 *   								Landcover (accessed from indiv) public members:
 *   - acflux_wood_harvest			gridcell-level C flux from harvest associated with cloning at landcover change (kgC/m2)
 *   - acflux_wood_harvest_lc[lc]	landcover-level C flux from harvest associated with cloning at landcover change (kgC/m2)
 *   - anflux_wood_harvest			gridcell-level N flux from harvest associated with cloning at landcover change (kgN/m2)
 *   - anflux_wood_harvest_lc[lc]	landcover-level N flux from harvest associated with cloning at landcover change (kgN/m2)
 */

void harvest_wood(Individual& indiv, double frac_cut, double harv_eff, double res_outtake_twig, 
		double res_outtake_coarse_root, bool lc_change) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	harvest_wood(indiv_cp, indiv.diam, indiv.pft, indiv.alive, frac_cut, harv_eff, res_outtake_twig, res_outtake_coarse_root);

	indiv_cp.copy_to_indiv(indiv, false, lc_change);

	if (lc_change) {
		Stand& stand = indiv.vegetation.patch.stand;
		Landcover& lc = stand.get_gridcell().landcover;
		lc.acflux_wood_harvest += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
		lc.anflux_wood_harvest += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
		if(stand.lc_origin < NLANDCOVERTYPES) {
			lc.acflux_wood_harvest_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
			lc.anflux_wood_harvest_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
		}
	}
}

/// Applies diameter rules if mt.diam_cut_low is set (default 0) and returns maximum man_strength for this individual. If not set, 1 is returned.
/** INPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following public members:
 *   - diam							tree diameter (m)
 *									ManagementType (accessed from indiv) public members:
 *	 - diam_cut_high				minimum tree diameter limit (cm) for cutting 100% of trees
 *   - firstmanageyear				calendar year when management starts (including suppression of disturbance and fire)
 *   - firstcutyear					first calendar year when wood harvest starts
 *   - firstcutyear_is_referenceyear	whether first cut year rather than patch age is used as a reference for timing of 
 *										cutting events
 *   - secondintervalstart			number of years after stand creation when the second (continuous) cutting interval starts
 *									Gridcellst (accessed from patch) public members:
 *   - diam_cut_low					minimum tree diameter limit (cm) for cutting (thinstrength*100)% of trees (can adapt to
 *									productivity)
 */
double diameter_rules(Individual& indiv) {

	// See Lagergren and Jönsson (2017) for the influence of site quality class of Swedish 
	// forests on diameter limits in target diameter cutting.

	Patch& patch = indiv.vegetation.patch;
	Gridcellst& gst = patch.stand.get_gridcell().st[patch.stand.stid];
	ManagementType& mt = patch.stand.get_current_management();

	double diam_cut_low = gst.diam_cut_low;	// mt diam_cut_low may be adjusted, stored in gst.diam_cut_low
	double diam_cut_high = mt.diam_cut_high;

	// Only use diameter limit for second (continuous cover) cutting period
	int first_manageyear = (mt.firstmanageyear < FAR_FUTURE_YEAR) ? mt.firstmanageyear - date.first_calendar_year : nyear_spinup;
	int first_cutyear = (mt.firstcutyear < FAR_FUTURE_YEAR) ? mt.firstcutyear - date.first_calendar_year : first_manageyear;
	int cutting_reference_age = mt.firstcutyear_is_referenceyear ? date.year - first_cutyear : patch.age;

	if(!diam_cut_low || mt.secondintervalstart == -1 || cutting_reference_age < mt.secondintervalstart)
		return 1.0;

	double man_strength;
	// indiv.diam is in m, diam_cut_low and diam_cut_high are in cm
	if (indiv.diam * 100 > diam_cut_low) {
		man_strength = patch.man_strength;
		if(indiv.diam * 100 > diam_cut_high)
			man_strength = 1.0;
	}
	else {
		man_strength = 0.0;
	}

	return man_strength;
}

/// Distributes man_strength (cutting intensity) from patch level to individual level
/** Calls check_harvest_cmass() to determine biomass at Individual and Patch levels.
 *  The patch management strength * cmass_wood demand is distributed to the individuals according to rules
 *  in this function and options in the function parameter list. The amount is conserved unless diameter_rules() returns a
 *  maximum man_strength value (individuals with diameters below limit left) and not enough wood cmass is available.
 *  In this case, the diameter limit is lowered by 1% each year until demand fulfilled.
 *
 *  INPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - man_strength 				management strength (cutting intensity)
 *  \param select_diam				Whether small (SELECT_DIAM_SMALL, 1) or large (SELECT_DIAM_LARGE, 2) diameter individuals
 *									are preferentially cut, trees above diam_limit only (SELECT_DIAM_LIMIT, 3).or no
 *									preference (SELECT_DIAM_NOPREF, 0); overrides select_age setting
 *  \param select_age				Whether young (SELECT_AGE_YOUNG, 1) or old (SELECT_AGE_OLD, 2) individuals are
 *									preferentially cut, or no preference (SELECT_AGE_NOPREF, 0)
 *  \param select_pft				Whether non-selected (SELECT_PFT_UNSEL, 1) or selected (SELECT_PFT_SEL, 2) pft:s are
 *									preferentially cut, unselected and selected cutting strengths specified separately
 *									(SELECT_PFT_UNSEL_SEL_SEPARATE, 3), shrubs and shade-intolerant pft:s preferentially cut
 *									(SELECT_PFT_SHADEINTOL, 4) or no preference (SELECT_PFT_NOPREF, 0)
 *  \param str_unsel				Cutting strength for unselected pft:s if select_pft = SELECT_PFT_UNSEL_SEL_SEPARATE
 *  \param str_sel					String of cutting strengths for selected pft:s if select_pft = SELECT_PFT_UNSEL_SEL_SEPARATE
 *									ManagementType (accessed from patch) public members:
 *   - firstmanageyear				calendar year when management starts (including suppression of disturbance and fire)
 *   - firstcutyear					first calendar year when wood harvest starts
 *   - firstcutyear_is_referenceyear	whether first cut year rather than patch age is used as a reference for timing
 *										of cutting events
 *   - secondcutinterval			wood cutting interval in years in the contiuous cutting period
 *   - secondintervalstart			number of years after stand creation when the second (continuous) cutting interval starts
 *   - adapt_diam_limit				whether to adapt diam_cut_low to forest stands with small trees
 *   - diam_cut_low					minimum tree diameter limit (cm) for cutting (thinstrength*100)% of trees
 *  INPUT/OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - distributed_cutting			whether this function has already been executed this year
 *									Gridcellst (accessed from patch) public members:
 *   - diam_cut_low					minimum tree diameter limit (cm) for cutting (thinstrength*100)% of trees (can adapt to productivity)
 *  OUTPUT PARAMETERS
 *									Individual (accessed from looping through patch vegetation) public members:
 *   - man_strength 				management strength (cutting intensity)
 */
void distribute_cutting(Patch& patch, int select_diam = SELECT_DIAM_NOPREF, int select_age = SELECT_AGE_NOPREF, 
	int select_pft = SELECT_PFT_NOPREF, double str_unsel = 0.0, double* str_sel = NULL) {

	/* Prevent individual man_strength to be set in several calls to this function the same year (called in manage_forest()
	 * before in set_forest_pft_structure()).
	 */
	if(patch.distributed_cutting)
		return;

	const int MAXAGE = 1500;
	Rank_individuals indiv_class(patch);
	ManagementType& mt = patch.stand.get_current_management();
	Stand& stand = patch.stand;
	StandType& st = stlist[stand.stid];
	const bool stem_cmass_only = true;
	double cmass_harvest_patch = check_harvest_cmass(patch, stem_cmass_only);
	double cmass_harvest_patch_selection = check_harvest_cmass(patch, stem_cmass_only, true);
	double cmass_harvest_remain = cmass_harvest_patch * patch.man_strength;
	double cmass_harvest_remain_init = cmass_harvest_remain;
	double cmass_harvest_patch_unsel = cmass_harvest_patch - cmass_harvest_patch_selection;

	double* cmass_pft = new double[stand.npft_selection];
	double* cmass_harvest_remain_pft = new double[stand.npft_selection];

	for(int i = 0; i < stand.npft_selection; i++) {
		cmass_pft[i] = 0.0;
		cmass_harvest_remain_pft[i] = 0.0;
	}

	// We need to re-calculate cmass_harvest_patch_selection for select_species case without some PFTs
	if(select_pft == SELECT_PFT_SHADEINTOL)
		cmass_harvest_patch_selection = 0.0;

	// Determine cmass:harvest for each cohort first for select_age != SELECT_AGE_NOPREF options
	double cmass_harvest_ageclass[MAXAGE + 1] = {0.0};
	double cmass_harvest_ageclass_selection[MAXAGE + 1] = {0.0};
	for(unsigned int i = 0; i < patch.vegetation.nobj; i++) {
		Individual& indiv = patch.vegetation[i];
		Standpft& spft = stand.pft[indiv.pft.id];
		if(indiv.age > MAXAGE)
			fail("Tree %d years old, increase MAXAGE in distribute_cutting()\n", indiv.age);
		bool pft_selection = spft.plant || mt.planting_system == "";
		if(select_pft == SELECT_PFT_SHADEINTOL) {
			// Exclude shrubs and shade-intolerant pfts from selection in this function
			pft_selection = pft_selection && !indiv.pft.is_shrub() && !indiv.pft.is_shade_intolerant_tree();
		}
		cmass_harvest_ageclass[(int)indiv.age] += check_harvest_cmass(indiv, stem_cmass_only);
		if(pft_selection) {
			double harvest_cmass_indiv = check_harvest_cmass(indiv, stem_cmass_only);
			if(select_pft == SELECT_PFT_SHADEINTOL)
				cmass_harvest_patch_selection += harvest_cmass_indiv;
			cmass_harvest_ageclass_selection[(int)indiv.age] += harvest_cmass_indiv;
			if(str_sel) {
				cmass_pft[stand.pft[indiv.pft.id].selection] += harvest_cmass_indiv;
				cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] += 
					str_sel[stand.pft[indiv.pft.id].selection] * harvest_cmass_indiv;
			}
		}
	}

	/* Two loops through the individuals of the patch if there is a preference in cutting selected or unselected pfts first;
	 * first lap for the preferred pfts, the second lap for the other pfts.
	 */
	int nlaps = (select_pft != SELECT_PFT_NOPREF) ? 2 : 1;

	for(int n=0; n<nlaps; n++) {

		double cmass_harvest_remain_init_selection;
		double cmass_harvest_remain_init_ageclass;
		int age_save = -1;	// Age of last individual
		bool cut_selection;

		if(select_pft == SELECT_PFT_UNSEL || select_pft == SELECT_PFT_SHADEINTOL) {
			if(!n) {
				cut_selection = false;
			}
			else {
				cut_selection = true;
			}
		}
		else {
			if(!n) {
				cut_selection = true;
			}
			else {
				cut_selection = false;
			}

			if(select_pft == SELECT_PFT_UNSEL_SEL_SEPARATE) {
				if(!n) { // Value used for selection in this case
					cmass_harvest_remain = cmass_harvest_patch_selection * patch.man_strength;
				}
				else {
					cmass_harvest_remain = cmass_harvest_patch_unsel * str_unsel;
				}
			}
		}

		for(unsigned int i = 0; i < patch.vegetation.nobj; i++) {

			if(cmass_harvest_remain < 1e-15 && !(cut_selection && select_pft == SELECT_PFT_UNSEL_SEL_SEPARATE && str_sel))
				break;

			int index;

			if(select_age == SELECT_AGE_YOUNG) {
				index = patch.vegetation.nobj - 1 - i;
			}
			else {
				index = i;
			}

			// select_diam overrides select_age
			if(select_diam == SELECT_DIAM_SMALL) {
				index = indiv_class.get_index(i);
			}
			else if(select_diam == SELECT_DIAM_LARGE) {
				index = indiv_class.get_index(patch.vegetation.nobj -1 - i);
			}

			if(!i)
				cmass_harvest_remain_init_selection = cmass_harvest_remain;

			Individual& indiv = patch.vegetation[index];

			if(indiv.pft.lifeform != TREE)
				continue;

			Standpft& spft = stand.pft[indiv.pft.id];
			bool pft_selection = spft.plant || mt.planting_system == "";

			if(select_pft == SELECT_PFT_SHADEINTOL) {
				// Exclude shrubs and shade-intolerant pfts from selection in this function
				pft_selection = pft_selection && !indiv.pft.is_shrub() && !indiv.pft.is_shade_intolerant_tree();
			}

			if(select_pft != SELECT_PFT_NOPREF) {
				if(!cut_selection && pft_selection)
					continue;
				if(cut_selection && !pft_selection)
					continue;
			}

			if(!i || (int)indiv.age != age_save)
				cmass_harvest_remain_init_ageclass = cmass_harvest_remain;
			age_save = (int)indiv.age;

			// Harvestable cmass_wood for this individual/cohort
			double cmass_harvest_cohort = check_harvest_cmass(indiv, stem_cmass_only);

			if(!cmass_harvest_cohort)
				continue;

			// select_diam overrides select_age
			if(select_diam != SELECT_DIAM_NOPREF) {

				// Determine if diameter rules restricts harvestable amount
				double max_cut = diameter_rules(indiv);

				// Satisfy cutting demand by cutting down each according to the preferences set by select_diam
				if(select_diam == SELECT_DIAM_LIMIT) {
					indiv.man_strength = max_cut;
				}
				else {
					if(cut_selection && select_pft == SELECT_PFT_UNSEL_SEL_SEPARATE && str_sel) {
						indiv.man_strength = min(1.0, cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] / cmass_harvest_cohort);
						cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] -= max(0.0, indiv.man_strength * cmass_harvest_cohort);
					}
					else {
						indiv.man_strength = min(1.0, cmass_harvest_remain / cmass_harvest_cohort);
						// if diam_cut_low set, don't cut smaller trees
						if(!max_cut)
							indiv.man_strength = 0.0;
					}
				}
			}
			else if(select_age != SELECT_AGE_NOPREF) {

				// Using equal cutting within an age-class
				if(select_pft != SELECT_PFT_NOPREF) {
					if(!cut_selection && (cmass_harvest_ageclass[(int)indiv.age] - cmass_harvest_ageclass_selection[(int)indiv.age])) {
						indiv.man_strength = min(1.0, cmass_harvest_remain_init_ageclass
											/ (cmass_harvest_ageclass[(int)indiv.age] - cmass_harvest_ageclass_selection[(int)indiv.age]));
					}
					else if(cut_selection && cmass_harvest_ageclass_selection[(int)indiv.age]) {
						if(select_pft == SELECT_PFT_UNSEL_SEL_SEPARATE && str_sel) {
							// ony one individual per pft per age !
							indiv.man_strength = min(1.0, cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] / cmass_harvest_cohort);
							cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] -= max(0.0, indiv.man_strength * cmass_harvest_cohort);
						}
						else {
							indiv.man_strength = min(1.0, cmass_harvest_remain_init_ageclass / cmass_harvest_ageclass_selection[(int)indiv.age]);
						}
					}
					else {
						indiv.man_strength = 0.0;
					}
				}
				else {
					if(cmass_harvest_ageclass[(int)indiv.age]) {
						indiv.man_strength = min(1.0, cmass_harvest_remain_init_ageclass / cmass_harvest_ageclass[(int)indiv.age]);
					}
					else {
						indiv.man_strength = 0.0;
					}
				}
			}
			else {

				if(select_pft != SELECT_PFT_NOPREF) {

					// First cut unwanted species by an equal amount
					if(!cut_selection && (cmass_harvest_patch - cmass_harvest_patch_selection)) {
						indiv.man_strength = min(1.0, cmass_harvest_remain_init_selection / (cmass_harvest_patch - cmass_harvest_patch_selection));
					}
					// then cut an equal amount of the pft:s in selection
					else if(cut_selection && cmass_harvest_patch_selection) {
						if(select_pft == SELECT_PFT_UNSEL_SEL_SEPARATE && str_sel) {
							indiv.man_strength = str_sel[stand.pft[indiv.pft.id].selection];	// Special case: not equal amounts !
						}
						else {
							indiv.man_strength = min(1.0, cmass_harvest_remain_init_selection / cmass_harvest_patch_selection);
						}
					}
					else {
						indiv.man_strength = 0.0;
					}
				}
				else {
					// Cut an equal amount of each cohort
					indiv.man_strength = patch.man_strength;
				}
			}
			cmass_harvest_remain -= indiv.man_strength * cmass_harvest_cohort;
		}
	}
	/* If diameter rules used, prescribed cutting may not be achieved (especially in young stands).
	 * Try to solve demand by reducing diameter limitby 1% each year when this happens.
	 */
	if(mt.secondcutinterval
		&& ((select_diam == SELECT_DIAM_LIMIT && cmass_harvest_remain_init && (cmass_harvest_remain_init - cmass_harvest_remain) < 1e-15) 
			|| (select_diam == SELECT_DIAM_SMALL || select_diam == SELECT_DIAM_LARGE) && cmass_harvest_remain > 1e-15)) {
		
		dprintf("Year %d: Warning: cmass_harvest_remain = %f, inital cmass_harvest_remain = %f; age = %d\n", 
			date.get_calendar_year(), cmass_harvest_remain,cmass_harvest_remain_init, patch.age);

		// Only in second (continuous) period
		int first_manageyear = (mt.firstmanageyear < FAR_FUTURE_YEAR) ? mt.firstmanageyear - date.first_calendar_year : nyear_spinup;
		int first_cutyear = (mt.firstcutyear < FAR_FUTURE_YEAR) ? mt.firstcutyear - date.first_calendar_year : first_manageyear;
		int cutting_reference_age = mt.firstcutyear_is_referenceyear ? date.year - first_cutyear : patch.age;

		if(mt.adapt_diam_limit && mt.diam_cut_low && mt.secondintervalstart > -1 && cutting_reference_age >= mt.secondintervalstart) {

			// Take into account number of patches that are cut each year:
			double harvests_per_year = (1.0 * patch.stand.npatch()) / mt.secondcutinterval;

			Gridcellst& gst = patch.stand.get_gridcell().st[patch.stand.stid];
			gst.diam_cut_low *= (1.0 - (0.01 / harvests_per_year));
			dprintf("New diam_cut_low = %f\n", gst.diam_cut_low);
		}
	}

	patch.distributed_cutting = true;

	if(cmass_pft)
		delete[] cmass_pft;
	if(cmass_harvest_remain_pft)
		delete[] cmass_harvest_remain_pft;

	return;
}

/// Sets management strength for individual trees to achieve prescribed tree pft composition
/** Calculates required fraction of tree pfts to cut to reach relative biomass target fractions specified in mt.targetfrac.
 *  Calls check_harvest_cmass() to determine pft and total biomass at Individual, Patch and Stand levels.
 *  Calls distribute_cutting() to set man_strength for tree individuals.
 *
 *  INPUT PARAMETERS
 *									ManagementType (accessed from looping through gridcell stand list) public members:
 *   - firsttargetyear				when to start cutting to reach pft target fractions (calendar year)
 *   - lasttargetyear				when to stop cutting to reach pft target fractions (calendar year)
 *   - targetfrac					string of pft target cmass fractions
 *   - targetcutinterval			interval of target cuttings
 *   - planting_system				rules for which pfts are planted
 *									Stand (accessed from looping through gridcell stand list) public members:
 *   - npft_selection				nNumber of pfts in selection
 *									Standpft (accessed from looping through pftlist) public members:
 *   - targetfrac					pft target cmass fraction for this StandType
 *   - selection					order of pft in selection string
 *									Patch (accessed from looping through stands) public members:
 *   - distributed_cutting			whether this function has already been executed this year
 *  OUTPUT PARAMETERS
 *									Individual (accessed from looping through patches' vegetation) public members:
 *   - man_strength 				management strength (cutting intensity)
 */
void set_forest_pft_structure(Gridcell& gridcell) {

	// Limit of sum of deviations of pft biomass fractions from target fractions that triggers cutting.
	const double DEVLIMIT = 0.1;

	for(unsigned int s=0;s<gridcell.nbr_stands();s++) {
		Stand& stand = gridcell[s];
		ManagementType& mt = stand.get_current_management();
		const bool stem_cmass_only = false;						// The two options will give small differences with different  
																// harvest_slow_frac (caused by rounding errors)
		const bool relax_target_cutting_after_pft_gone = false; // Default value false will prioritise pft selection over 
																// stand age (all other trees may be cut when one pft dies)
		const bool no_target_cutting_before_man_start = true;	// No target-cutting before management starts

		bool management_started = stand[0].managed;				// patch.managed is the same for all patches in a stand

		int first_targetyear = nyear_spinup; // Simulation year when target cutting starts; default is directly after spinup.
		if(mt.firsttargetyear < FAR_FUTURE_YEAR)	// Initialised to FAR_FUTURE_YEAR; other values set in instruction file.
			first_targetyear = mt.firsttargetyear - date.first_calendar_year;
		if(mt.planting_system != "SELECTION" || (mt.targetfrac == "" && !readtargetcutting)
			|| (readtargetcutting && mt.targetfrac_input_mode == 2)
			|| date.get_calendar_year() > mt.lasttargetyear || date.year < first_targetyear
			|| (no_target_cutting_before_man_start && !management_started)) {

			continue;
		}

		double* target = new double[stand.npft_selection];
		for(int i = 0; i < stand.npft_selection; i++) {
			target[i] = 0.0;
		}
		double target_sum = 0.0;

		pftlist.firstobj();
		while (pftlist.isobj) {
			Pft& pft = pftlist.getobj();
			Standpft& spft = stand.pft[pft.id];
			if(mt.pftinselection((const char*)pft.name) && spft.selection != -1) {
				target[spft.selection] = spft.targetfrac;
				target_sum += spft.targetfrac;
			}
			pftlist.nextobj();
		}

		// When target sum is zero for this year (in input file), no target cutting will occur.
		if(!target_sum)
			continue;

		// Normalise targets if target sum > 1.0:
		for(int i=0;i<stand.npft_selection;i++) {
			if(target_sum > 1.0)
				target[i] /= target_sum;
		}

		// Check deviations at stand level:
		double target_sum_stand = 0.0;
		double* target_stand = new double[stand.npft_selection];
		double* cmass_pft_stand = new double[stand.npft_selection];
		double* cutstr_pft_stand = new double[stand.npft_selection];

		for(int i = 0; i < stand.npft_selection; i++) {
			target_stand[i] = 0.0;
			cmass_pft_stand[i] = 0.0;
			cutstr_pft_stand[i] = 0.0;
		}

		double cmass_unselected_stand = 0.0;
		double cmass_total_stand = check_harvest_cmass(stand, stem_cmass_only);

		if(!cmass_total_stand)
			continue;

		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			for(unsigned int i = 0; i < patch.vegetation.nobj; i++) {
				Individual& indiv = patch.vegetation[i];
				Standpft& spft = stand.pft[indiv.pft.id];
				if(mt.pftinselection((const char*)indiv.pft.name) && spft.selection != -1) {
					cmass_pft_stand[spft.selection] += check_harvest_cmass(indiv, stem_cmass_only) / stand.nobj;
				}
				else if(indiv.pft.lifeform == TREE) {
					cmass_unselected_stand += check_harvest_cmass(indiv, stem_cmass_only) / stand.nobj;
				}
			}
			stand.nextobj();
		}

		// Avoid excessive cutting after species disappearance
		double exclude_frac_stand = 0.0;
		pftlist.firstobj();
		while(pftlist.isobj) {
			Pft& pft = pftlist.getobj();
			Standpft& spft = stand.pft[pft.id];
			if(spft.selection != -1 && target[spft.selection] > 0.0) {
				if(!cmass_pft_stand[spft.selection] && relax_target_cutting_after_pft_gone) {
					exclude_frac_stand += target[spft.selection];
				}
				else {
					target_stand[spft.selection] = target[spft.selection];
				}
			}
			pftlist.nextobj();
		}
		if(exclude_frac_stand) {
			for(int i=0;i<stand.npft_selection;i++) {
				if(1.0 - exclude_frac_stand)
					target_stand[i] /= (1.0 - exclude_frac_stand);
			}
		}

		double remove_cmass_selected_stand = 0.0;
		double cutstr_selected_stand = 0.0;

		for(int i=0;i<stand.npft_selection;i++) {
			double remove_cmass_stand = 0.0; 
			if(target_stand[i] - 1.0) {
				if((cmass_pft_stand[i] / cmass_total_stand) > 0.99 && target_stand[i] <= 0.99 && relax_target_cutting_after_pft_gone) {
					dprintf("Target cutting suspended year %d: cmass_pft_stand[i] / cmass_total_stand = %f\n\n",
						date.get_calendar_year(), cmass_pft_stand[i] / cmass_total_stand);
				}
				else {
					remove_cmass_stand = max(0.0, (target_stand[i] * cmass_total_stand - cmass_pft_stand[i]) / (target_stand[i] - 1.0));
				}
			}

			remove_cmass_selected_stand += remove_cmass_stand;
			if(cmass_pft_stand[i])
				cutstr_pft_stand[i] += remove_cmass_stand / cmass_pft_stand[i];
			cutstr_selected_stand += cutstr_pft_stand[i];
			target_sum_stand += target_stand[i];
		}

		double fraction_unselected_stand = cmass_unselected_stand / cmass_total_stand;
		double cutstr_unselected_stand = 0.0;
		double remove_cmass_unselected_stand = 0.0;
		if(target_sum_stand)
			remove_cmass_unselected_stand = max(0.0, ((1.0 - target_sum_stand) * cmass_total_stand - cmass_unselected_stand) / -target_sum_stand);
		if(cmass_unselected_stand)
			cutstr_unselected_stand = remove_cmass_unselected_stand / cmass_unselected_stand;

		// Equal importance of pft relative deviations in- and outside of selection
		double cutstr_total_stand = cutstr_selected_stand + cutstr_unselected_stand;
		// Importance of pft relative deviations in- and outside of selection weighted by cmass of selected and unselected pfts.
		// cutstr_total_stand = (remove_cmass_selected_stand + remove_cmass_unselected_stand) / cmass_total_stand;

		// Check deviations at patch level and set man_strength>:
		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();

			// Suppress target cutting during secondary period if mt.suppress_second_target = true
			bool suppress_secondary = false;
			int first_manageyear = (mt.firstmanageyear < FAR_FUTURE_YEAR) ? mt.firstmanageyear - date.first_calendar_year : nyear_spinup;
			int first_cutyear = (mt.firstcutyear < FAR_FUTURE_YEAR) ? mt.firstcutyear - date.first_calendar_year : first_manageyear;
			int cutting_reference_age = mt.firstcutyear_is_referenceyear ? date.year - first_cutyear : patch.age;
			if(mt.suppress_second_target && mt.secondintervalstart > -1 && cutting_reference_age >= mt.secondintervalstart)
				suppress_secondary = true;

			if(!(patch.age >= mt.targetstartage && !(patch.age % mt.targetcutinterval) && !patch.distributed_cutting
					&& !patch.clearcut_this_year && !suppress_secondary)) {
				stand.nextobj();
				continue;
			}

			double target_sum_patch = 0.0;

			double* cmass_pft = new double[stand.npft_selection];
			double* cutstr_pft = new double[stand.npft_selection];
			double* target_patch = new double[stand.npft_selection];

			for(int i = 0; i < stand.npft_selection; i++) {
				cmass_pft[i] = 0.0;
				cutstr_pft[i] = 0.0;
				target_patch[i] = 0.0;
			}

			double cmass_unselected = 0.0;
			double cmass_total = check_harvest_cmass(patch, stem_cmass_only);

			if(!cmass_total) {
				stand.nextobj();
				continue;
			}

			for(unsigned int i = 0; i < patch.vegetation.nobj; i++) {
				Individual& indiv = patch.vegetation[i];
				Standpft& spft = stand.pft[indiv.pft.id];
				if(mt.pftinselection((const char*)indiv.pft.name) && spft.selection != -1) {
					cmass_pft[spft.selection] += check_harvest_cmass(indiv, stem_cmass_only);
				}
				else if(indiv.pft.lifeform == TREE) {
					cmass_unselected += check_harvest_cmass(indiv, stem_cmass_only);
				}
			}

			// Avoid ezcessive cutting after species disappearance
			double exclude_frac_patch = 0.0;
			pftlist.firstobj();
			// Loop through PFTs
			while(pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				Standpft& spft = stand.pft[pft.id];
				if(spft.selection != -1 && target[spft.selection] > 0.0) {		
					if(!cmass_pft[spft.selection] && relax_target_cutting_after_pft_gone) {
						exclude_frac_patch += target[spft.selection];
					}
					else {
						target_patch[spft.selection] = target[spft.selection];
					}
				}
				pftlist.nextobj();
			}
			if(exclude_frac_patch) {
				for(int i=0;i<stand.npft_selection;i++) {
					if(1.0 - exclude_frac_patch)
						target_patch[i] /= (1.0 - exclude_frac_patch);
				}
			}

			double remove_cmass_selected = 0.0;
			double cutstr_selected = 0.0;

			for(int i=0;i<stand.npft_selection;i++) {
				double remove_cmass = 0.0;
 
				if(target_patch[i] - 1.0) {
					if((cmass_pft[i] / cmass_total) > 0.99 && target_patch[i] <= 0.99 && relax_target_cutting_after_pft_gone) {
						dprintf("Target cutting suspended year %d: cmass_pft[i] / cmass_total = %f\n\n",
							date.get_calendar_year(), cmass_pft[i] / cmass_total);
					}
					else {
						remove_cmass = max(0.0, (target_patch[i] * cmass_total - cmass_pft[i]) / (target_patch[i] - 1.0));
					}
				}

				remove_cmass_selected += remove_cmass;
				if(cmass_pft[i])
					cutstr_pft[i] = remove_cmass / cmass_pft[i];
				cutstr_selected += cutstr_pft[i];
				target_sum_patch += target_patch[i];
			}

			double fraction_unselected = cmass_unselected / cmass_total;
			double cutstr_unselected = 0.0;
			double remove_cmass_unselected = 0.0;
			if(target_sum_patch)
				remove_cmass_unselected = max(0.0, ((1.0 - target_sum_patch) * cmass_total - cmass_unselected) / -target_sum_patch);
			if(cmass_unselected)
				cutstr_unselected = remove_cmass_unselected / cmass_unselected;

			// Equal importance of pft relative deviations in- and outside of selection
			double cutstr_total = cutstr_selected + cutstr_unselected;
			// Importance of pft relative deviations in- and outside of selection weighted by cmass of selected and unselected pfts.
			// cutstr_total = (remove_cmass_selected + remove_cmass_unselected) / cmass_total;

			/*  Modes of cutting:
			 * 1: Cut when patch fraction deviations > DEVLIMIT, use patch overshoot values (default)
			 * 2: Cut when stand fraction deviations > DEVLIMIT, use patch overshoot values 
			 * 3: Cut when stand fraction deviations > DEVLIMIT, use stand overshoot values
			 */
			int targetcutmode = mt.targetcutmode;

			double cutstr_total_use = cutstr_total;
			double cutstr_unselected_use = cutstr_unselected;
			double *cutstr_pft_use = cutstr_pft;

			if(targetcutmode > 1)
				cutstr_total_use = cutstr_total_stand;
			if(targetcutmode == 3) {
				cutstr_unselected_use = cutstr_unselected_stand;
				cutstr_pft_use = cutstr_pft_stand;
			}
	
			if(cutstr_total_use > DEVLIMIT)
				distribute_cutting(patch, mt.targetthinselectdiam, mt.targetthinselectage, SELECT_PFT_UNSEL_SEL_SEPARATE,
									cutstr_unselected_use, cutstr_pft_use);	// Default: cut largest trees first

			if(cmass_pft)
				delete[] cmass_pft;
			if(cutstr_pft)
				delete[] cutstr_pft;
			if(target_patch)
				delete[] target_patch;

			stand.nextobj();
		}
		if(target)
			delete[] target;
		if(target_stand)
			delete[] target_stand;
		if(cmass_pft_stand)
			delete[] cmass_pft_stand;
		if(cutstr_pft_stand)
			delete[] cutstr_pft_stand;
	}
}

/// Returns true if tree density is below clear-cutting limit mt.dens_target_cc when ifclearcut_by_density = true.
/** INPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - age							patch age
 *									ManagementType (accessed from patch) public members:
 *   - planting_system				rules for which pfts are planted; custom pft group planting_systems are used in
 *									this function (optionally)
 *   - firstclearcutyear			first calendar year when clear-cut allowed
 *   - delayduecutting				number of years to distribute clearcut of patches that were due to be cut before
 *									firstclearcutyear
 *   - dens_target_cc				stand density lower limit that triggers clearcut
 *  INPUT/OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - cut_due						whether tree density is below limit at firstclearcutyear, reset clear-cut year
 */
bool clearcut_by_density(Patch& patch) {

	// Based on clearcut scheme in Bellassen (2010), but tree density and density parameters were selected to give 
	// rotation times around 100 years in the early 2000s in LPJ-GUESS simulations.

	ManagementType& mt = patch.stand.get_current_management();
	StandType& st = stlist[patch.stand.stid];

	if(date.get_calendar_year() < mt.firstclearcutyear - 1)
		return false;

	double dens_target_cc = 150.0;	// Default value if not set below

	// Default values for tree functional groups, can be changed by mt values below
	if(mt.planting_system == "NEEDLELEAF_EVERGREEN" || mt.planting_system == "NEEDLELEAF_DECIDUOUS") {
		dens_target_cc = DENSTARGET_NL;
	}
	else if(mt.planting_system == "BROADLEAF_EVERGREEN" || mt.planting_system == "BROADLEAF_DECIDUOUS") {
		dens_target_cc = DENSTARGET_BL;
	}

	if(mt.dens_target_cc)
		dens_target_cc = mt.dens_target_cc;

	double dens = 0.0;
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (indiv.pft.lifeform == TREE) {
			dens += indiv.densindiv;
		}
		vegetation.nextobj();
	}

	// Convert from m-2 to ha-1
	dens *= M2_PER_HA;

	if(patch.age > 10 && dens < dens_target_cc && dens) {

			if(date.get_calendar_year() == mt.firstclearcutyear - 1)
				patch.cut_due = true;

			if(date.get_calendar_year() >= mt.firstclearcutyear &&
				(date.get_calendar_year() >= mt.firstclearcutyear + mt.delayduecutting
					|| !patch.cut_due || !(patch.age % mt.delayduecutting))) {

				patch.cut_due = false;
				return true;
			}
	}

	return false;
}

/// Setting of initial density when using thin_reineke (to avoid dependence on first_cutyear value)
/** INPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - age							patch age; reset at clear-cut
 *									Individual (accessed from looping through patches' vegetation) public members:
 *   - densindiv 					average density of individuals over patch (indiv/m2)
 *  OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - dens_start					initial tree individual density
 */
void thin_reineke_init(Patch& patch) {

	if(patch.age != 1)
		return;

	double dens = 0.0;

	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		if (indiv.pft.lifeform == TREE)
			dens += indiv.densindiv;
		vegetation.nextobj();
	}
	// Convert from m-2 to ha-1
	dens *= M2_PER_HA;

	patch.dens_start = dens;
}

/// Automated thinning routine based on Reineke's self-thinning rule
/** Self-thinning is avoided by cutting when the quadratic mean diameter (diam_sq) approaches the self-thinning curve,
 *  which is derived from diam_sq-density log-log plots of model output.
 *  Sets patch.man_strength and calls distribute_cutting() to set man_strength for individuals.
 *  
 *  INPUT PARAMETERS
 *  \param patch					reference to a Patch
 *									ManagementType (accessed from patch) public members:
 *   - planting_system				rules for which pfts are planted; custom pft group planting_systems are used in this
 *									function (optionally)
 *   - alpha_st						self-thinning parameter (theoretical density when Dg = 1m using the density-Dg equation)
 *   - rdi_target					thinning "intensity" (low value more intense) 
 *   - dens_target_cc				stand density lower limit that triggers clearcut
 *									Individual (accessed from looping through patch vegetation) public members:
 *   - densindiv 					average density of individuals over patch (indiv/m2)
 *   - diam		 					stem diameter (m)
 *  OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - harvest_to_litter			switch for non-commercial thinning (harvest goes to litter)
 */
void thin_reineke(Patch& patch) {

	/* Based on Reineke's self-thinning rule, as in Bellassen (2010)
	 *
	 * dens_max = alpha_st / pow(Dg, beta_st), where dens_max is stand maximum density before self-thinning (ind/ha),
	 * alpha_st and beta_st are parameters and Dg is the quadratic mean diameter (m),
	 * Dg = pow(diameter square sum, 0.5) (equation in Bellassen 2010 is incorrect)
	 * 
	 * alpha_st and beta_st were calibrated from log-log plots of quadratic mean diameter (Dg) and tree density:
	 * log Dg=log alpha_st/beta_st-1/beta_st*log dens and density parameters were selected to give rotation times
	 * around 100 years in the early 2000s in LPJ-GUESS simulations.
	 *
	 * To avoid self-thinning mortality, the relative density index (rdi) = dens / dens_max, is monitored
	 * (dens_max = alpha_st / pow(Dg, beta_st) and kept close to a target value, rdi_target by cutting when
	 * rdi reaches (rdi_target + delta_rdi) to reach (rdi_target - delta_rdi)
	 */

	ManagementType& mt = patch.stand.get_current_management();
	StandType& st = stlist[patch.stand.stid];

	double alpha_st = 0.0;		// Default value set below if not specified in instruction file management parameters
	double beta_st = 1.6;
	// thinning "intensity" (low value more intense)
	double rdi_target = 0.75;	// Bellassen value, modified below
	double dens_start = 5000.0;	// Changed this from Bellassen's 10000 to represent what we use in LPJ-GUESS if not 
								// specifically specifying a higher planting density
	double dens_target = 0.0;	// Default value set below if not specified in instruction file management parameters

	// Default values for tree functional groups, can be changed by mt values below
	if(mt.planting_system == "NEEDLELEAF_EVERGREEN" || mt.planting_system == "NEEDLELEAF_DECIDUOUS") {

		alpha_st = 65;					// Values obtained from LPJ-GUESS simulations of needle-leaf monocultures without
										// re-establishment (60: shorter rotation time)
		rdi_target = 0.7;				// Bellassen value 0.75; changed to 0.7 to get shorter rotation time (with alpha_st 65)
		dens_target = DENSTARGET_NL;
	}
	else if(mt.planting_system == "BROADLEAF_EVERGREEN" || mt.planting_system == "BROADLEAF_DECIDUOUS") {

		alpha_st = 40;					// Values obtained from LPJ-GUESS simulations of broad-leaf monocultures and mixed
										// stands without re-establishment
		rdi_target = 0.85;				// 0.85-0.9; 0.85 for less mortality
		dens_target = DENSTARGET_BL;
	}

	if(mt.alpha_st)
		alpha_st = mt.alpha_st;
	if(mt.rdi_target)
		rdi_target = mt.rdi_target;
	if(mt.dens_target_cc)
		dens_target = mt.dens_target_cc;

	if(!alpha_st) {
		alpha_st = 65;			// Using needle-leaf value if value not set at this stage.
		dprintf("Warning: alpha_st value not set for %s, using needle-leaf value\n", (char*)st.name);
	}
	if(!dens_target) {
		dens_target = 150.0;	// Using default value if value not set at this stage.
		dprintf("Warning: dens_target value not set for %s, using average default value\n", (char*)st.name);
	}

	double dens = 0.0;
	double diam_sq = 0.0;

	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (indiv.pft.lifeform == TREE) {
			dens += indiv.densindiv;
			// Weight the different cohorts
			diam_sq += indiv.diam * indiv.diam * indiv.densindiv;
		}
		vegetation.nextobj();
	}

	diam_sq /= dens;

	// Convert from m-2 to ha-1
	dens *= M2_PER_HA;

	double Dg = pow(diam_sq, 0.5);
	double dens_max = alpha_st / pow(Dg, beta_st);
	double rdi = dens / dens_max;

	if(patch.dens_start)
		dens_start = patch.dens_start;
	double delta_rdi = min(0.1, 0.05 + 0.05 * max(0.0, log(dens / dens_target) / log(dens_start / dens_target)));

	if(rdi > (rdi_target + delta_rdi)) {
		patch.man_strength = (rdi - (rdi_target - delta_rdi)) / rdi;
		// Pre-commercial thinning: harvested biomass to litter
		if(patch.age < 20)
			patch.harvest_to_litter = true;
		// Cut young trees first & cut shrubs and shade-intolerant species first.
		distribute_cutting(patch, SELECT_DIAM_NOPREF, SELECT_AGE_YOUNG, SELECT_PFT_SHADEINTOL);
	}
}

/// Performs forest management in all stands this year
/** Sets patch.man_strength in calls to manage_forest() and set_forest_pft_structure(). Harvest is done in calls to harvest_forest().
 *  Sets new managements in calls to forest_rotation().
 *  Methods described in Lindeskog et al. 2021.
 */
void manage_forests(Gridcell& gridcell) {

	if (!run_landcover || date.day) {
		return;
	}

	// Determine how much to harvest in all forest and natural stands.
	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;

		// Avoid calling manage_forest() when converting harvest biomass input to man_strength
		// (harvest_secondary_to_new_stand = false, not currently implemented).
		if(harvest_secondary_to_new_stand && (stand.landcover == FOREST || stand.landcover == NATURAL)) {
			
			stand.firstobj();
			while (stand.isobj) {
				Patch& patch = stand.getobj();
				if(harvest_secondary_to_new_stand)	{	
					manage_forest(patch);
				}
				stand.nextobj();
			}
		}

		forest_rotation(stand);

		++gc_itr;
	}

	// Target cutting; bypassed in years when normal cuttings are set in manage_forest()
	set_forest_pft_structure(gridcell);

	// Perform the previously determined harvest for all the individuals in all forest and natural stands.
	Gridcell::iterator gc_itr2 = gridcell.begin();
	while (gc_itr2 != gridcell.end()) {
		Stand& stand = *gc_itr2;

		stand.firstobj();
		while (stand.isobj && (stand.landcover == FOREST || stand.landcover == NATURAL)) {
			Patch& patch = stand.getobj();
			Vegetation& vegetation = patch.vegetation;
			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();

				bool killed = false;
				harvest_forest(indiv, indiv.pft, indiv.alive, 0.0, killed);

				if(!killed)
					vegetation.nextobj();
			}
			stand.nextobj();
		}
		++gc_itr2;
	}
}

/// Sets biomass fraction to be cut in clear-cuts and thinnings for all trees in this patch this year.
/** Sets patch.man_strength (tree biomass to be cut, 0-1) in clearcut- and continuous management schemes (harvest_system "CLEARCUT" and "CONTINUOUS").
 *  Calls distribute_cutting() to set man_strength for individuals in thinnings.
 *  Sets patch.man_strength to 1 in first_manageyear if mt.cutfirstyear == true in stands created at the start of the simulation.
 *  Sets patch.managed in first_manageyear (disables disturbance and fire if suppress_disturbance or suppress_fire == true).
 *  INPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - age			 				patch age
 *   - man_strength 				management strength (cutting intensity)
 *									ManagementType (accessed from patch) public members:
 *   - cutfirstyear					whether the patch should be clear-cut firstmanageyear
 *   - firstmanageyear				calendar year when management starts (including suppression of disturbance and fire)
 *   - firstcutyear					first calendar year when wood harvest starts
 *   - firstclearcutyear			first calendar year when clear-cut allowed
 *   - ifthin_reineke				whether automated thinning according to self-thin rule is selected
 *   - firstcutyear_is_referenceyear	whether first cut year rather than patch age is used as a reference for timing
 *										of cutting events
 *   - secondintervalstart			number of years after stand creation when the second (continuous) cutting interval starts
 *   - ifclearcut_by_density		whether clear-cut triggered by tree density going below a limit is selected
 *   - ifclearcut_optimal_age		whether clear-cut triggered by growth reaching diminishing return is selected
 *   - distribute_cuttings_among_patches		whether cuttings are evenly distributed among patches
 *   - harvest_system				whether clear-cut or continuous cutting schemes are selected
 *   - cutinterval					rotation time in years
 *   - thinstrength[][]				cutting strength of thinning events
 *   - thinstrength_unsel[][]		cutting strength for unselected pft:s if select_pft = 3
 *   - thintime[][]					timing of thinning events, relative to rotation period
 *   - thinselectdiam[][]			whether small (1) or large (2) diameter individuals are preferentially cut, trees above.
 *									diam_limit only (3).or no preference (0)
 *   - thinselectage[][]			whether young (1) or old (2) individuals are preferentially cut, or no preference (0)
 *   - thinselectpft[][]			whether non-selected (1) or selected (2) pft:s are preferentially cut, unselected and 
 *									selected cutting strengths specified separately (3), shrubs and shade-intolerant
 *									pft:s preferentially cut (4) or no preference (0)
 *									Stand (accessed from patch) public members:
 *   - first_year					simulation year when this stand was created
 *									Gridcellst (accessed from patch) public members:
 *   - cutinterval_st				cutting interval read from input file
 *									Global parameters:
 *   - readcutinterval_st			whether to use cutting interval in input file
 *  OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - age							patch age; reset at clear-cut
 *   - man_strength					patch-level management strength (cutting intensity)
 *   - managed						whether management has started (including suppression of disturbance and fire)
 *   - plant_this_year				whether patch should be planted this year
 *   - clearcut_this_year			whether patch has been clear-cut this year
 *   - has_been_cut					whether patch has ever been cut
 *   - cutinterval_actual			number of years since last clear-cut
 *   - harvest_to_litter			switch for non-commercial thinning (harvest goes to litter)
 *									Gridcellst (accessed from patch) public members:
 *   - reset_cutinterval_st			switch to reset cutinterval_st this year
 */
void manage_forest(Patch& patch) {

	Stand& stand = patch.stand;
	StandType& st = stlist[stand.stid];
	ManagementType& mt = stand.get_current_management();
	Gridcellst& gcst = stand.get_gridcell().st[st.id];

	int first_manageyear = nyear_spinup; // Simulation year when forestry management starts; default directly after spinup.

	if(mt.firstmanageyear < FAR_FUTURE_YEAR)	// Initialised to FAR_FUTURE_YEAR; other values set in instruction file.
		first_manageyear = mt.firstmanageyear - date.first_calendar_year;

	// firstmanageyear from input file overrides mt firstmanageyear value
	if(readfirstmanageyear_st && gcst.firstmanageyear_st < FAR_FUTURE_YEAR)
		first_manageyear = gcst.firstmanageyear_st - date.first_calendar_year;

	if(date.year < first_manageyear || !mt.is_managed())
		return;

	patch.managed = true;

	if(patch.stand.first_year == date.year) {
		patch.plant_this_year = true;	// Forces establishment first stand year to behave like after clearcut
		patch.has_been_cut = true;
	}

	double cut_fraction = 0.0;
	double cut_fraction_unsel = 0.0;
	int cut_interval = mt.cutinterval;
	bool clearcut_now = false;

	// clearcut_first_manageyear: clearcut in year==first_manageyear only in stands created at start of simulation (false),
	//							  clearcut in year==first_manageyear in all stands (true)
	bool clearcut_first_manageyear = false;
	if(mt.cutfirstyear && date.year == first_manageyear && (stand.first_year == 0 || clearcut_first_manageyear))
		clearcut_now = true;

	if(mt.ifthin_reineke)
		thin_reineke_init(patch);

	int first_cutyear = first_manageyear; // Simulation year when forestry harvesting starts; default is same year as 
										  // first_manageyear.
	if(mt.firstcutyear < FAR_FUTURE_YEAR) // Initialised to 1000000; other values set in instruction file.
		first_cutyear = mt.firstcutyear - date.first_calendar_year;

	/* If mt.firstcutyear_is_referenceyear = true, years since first_cutyear rather than patch age is used as a reference
	 * for timing of cutting events
	 */
	int cutting_reference_age = mt.firstcutyear_is_referenceyear ? date.year - first_cutyear : patch.age;

	// cutinterval from input file overwrites mt cutinterval value and other clearcut triggers if value for st exists in
	// the input file.
	if(readcutinterval_st && gcst.cutinterval_st != 0)
		cut_interval = (int)gcst.cutinterval_st;

	/* Cut according to number of patches and patch id to get an even patch age distribution.
	 * Number of patches (npatch_secondarystand) and patch id determines when patch is clearcut (age reset to 0) during
	 * a period of cutinterval (harvest_system CLEARCUT) or secondintervalstart (harvest_system CONTINUOUS) years after
	 * stand creation.
	 */
	bool in_distribute_patch_ages_period = false;
	int nyears_distribute_patch_ages = (mt.harvest_system == "CONTINUOUS") ? mt.secondintervalstart : cut_interval;

	if(mt.distribute_patch_ages && nyears_distribute_patch_ages
		&& date.year < max(stand.first_year, stand.clone_year) + nyears_distribute_patch_ages) {

		int patch_order = (int)(patch.id * nyears_distribute_patch_ages * 1.0 / (1.0 * stand.npatch()));

		/* patch age cutting overrides first (regrowth) period cutting in continuous harvest and thinnings in the first
		 * rotation period in clearcut.
		 */
		in_distribute_patch_ages_period = true;

		if(!((date.year - max(stand.first_year, stand.clone_year) - patch_order) % nyears_distribute_patch_ages)) {
			if(date.year - max(stand.first_year, stand.clone_year) > 0)
				clearcut_now = true;
		}
	}

	if(date.year < first_cutyear && !clearcut_now)
		return;

	if(mt.harvest_system == "CLEARCUT") {

		if(patch.age && !clearcut_now && !in_distribute_patch_ages_period) {

			/* Clearcut: priority when several cutting methods defined in the instruction file: 
			 * cut_interval input file > ifclearcut_by_density > cut_interval defined in st/mt > ifclearcut_optimal_age
			 */

			// Use density limit to trigger clearcut
			if(mt.ifclearcut_by_density && !readcutinterval_st) { 
				if(clearcut_by_density(patch)) {
					clearcut_now = true;
				}
			}
			// Clearcut rotation age defined in instruction file or input file
			else if(cut_interval) {
				// Distribute clearcuts among patches
				if(mt.distribute_cuttings_among_patches) {
					int patch_order = (int)(patch.id * cut_interval * 1.0 / (1.0 * stand.npatch()));
					if(!((date.year - max(stand.first_year, stand.clone_year) - patch_order) % cut_interval))
					// if(!((date.year - first_cutyear - patch_order) % cut_interval))	// patch 0 wil be cut firstcutyear
																						// (synchronised cuttings in all stands)
						clearcut_now = true;
				}
				else if(readcutinterval_st) {
					if(cutting_reference_age >= cut_interval) {
						clearcut_now = true;
					}
				}
				// Use original patch ages
				else if (!(cutting_reference_age % cut_interval)) {
					clearcut_now = true;
				}
			}
			// Use optimal rotation age to trigger clearcut
			else if(mt.ifclearcut_optimal_age) {
				/* First attempt to calculate optimum rotation age for clearcut.
				 * Clearcut occurs when last five year patch wood increment is dropping below overall achieved
				 * annual wood increment (excluding shrubs), indicating a diminishing return.
				 */
				if(patch.cmass_wood(true) / max(1, patch.age) > patch.get_tree_cmass_wood_inc_5() && patch.age > 20) {
					clearcut_now = true;
				}
			}

			if(date.get_calendar_year() < mt.firstclearcutyear) {
				clearcut_now = false;
			}

			// Thinnings

			if(!clearcut_now)  {
				// Thinnings according to Reineke's rule
				if(mt.ifthin_reineke) {
					thin_reineke(patch);
				}
				// Thinnings defined in instruction file
				else if(cut_interval) {
					int age = cutting_reference_age;
					// Distribute continuous cuttings among patches (following clear-cut distribution)
					if(mt.distribute_cuttings_among_patches) {
						int patch_order = (int)(patch.id * cut_interval * 1.0 / (1.0 * stand.npatch()));
						age = date.year - max(stand.first_year, stand.clone_year) - patch_order;
						// age = date.year - first_cutyear - patch_order;	// patch 0 wil be clear-cut firstcutyear
																			// (synchronised cuttings in all stands)
					}
					for(int t=0;t<NTHINNINGS;t++) {
						if((mt.thinstrength[0][t] || mt.thinstrength_unsel[0][t])
							&& (age == (int)(floor((cut_interval * mt.thintime[0][t]) + 0.5)))) {

							cut_fraction = mt.thinstrength[0][t];
							cut_fraction_unsel = mt.thinstrength_unsel[0][t];
							patch.man_strength = cut_fraction;
							// Pre-commercial thinning: harvested biomass to litter
							if(!t)
								patch.harvest_to_litter = true;
							distribute_cutting(patch, mt.thinselectdiam[0][t], mt.thinselectage[0][t], mt.thinselectpft[0][t], cut_fraction_unsel);
						}
					}
				}
			}
		}
	}

	if(clearcut_now) {
		cut_fraction = 1.0;
		patch.man_strength = cut_fraction;
		patch.cutinterval_actual = patch.age;
		patch.cutinterval_actual_thisyear = patch.age;
		patch.age = 0;
		patch.plant_this_year = true;
		patch.clearcut_this_year = true;
		gcst.reset_cutinterval_st = true;
	}
	else if(mt.harvest_system == "CONTINUOUS" && !in_distribute_patch_ages_period) {

		if(!cut_interval) {

			/* See Lagergren and Jönsson (2017) for the influence of site quality class (sqc) of Swedish 
			 * forests on rotation period in forest management.
			 */
			const double sqc_min = 2.351;	// The minimum average sqc for a county in Sweden
			const double sqc_max = 11.311;	// The maximum average sqc" for a county in Sweden
			const double sqc = 10.0;		// Temporary static value (gives a cut_interval of 17 years)

			cut_interval=30-(int)(15.0*(sqc-sqc_min)/(sqc_max-sqc_min));
		}

		int n = 0;	// thinningloop
		int age = cutting_reference_age;
		if(mt.secondintervalstart > -1 && cutting_reference_age >= mt.secondintervalstart) {
			n = 1;
			cut_interval = mt.secondcutinterval;
			age = cutting_reference_age - mt.secondintervalstart;

			// Distribute continuous cuttings among patches
			if(mt.distribute_cuttings_among_patches) {
				int patch_order = (int)(patch.id * cut_interval * 1.0 / (1.0 * stand.npatch()));
				age = date.year - max(stand.first_year, stand.clone_year) - patch_order;
				// age = date.year - first_cutyear - patch_order;	// patch 0 wil be cut firstcutyear
			}														// (synchronised cuttings in all stands)
		}

		for(int t=0;t<NTHINNINGS;t++) {

			if((mt.thinstrength[n][t] || mt.thinstrength_unsel[n][t])
				&& (age % cut_interval) == (int)(floor((cut_interval * mt.thintime[n][t]) + 0.5))) {

				cut_fraction = mt.thinstrength[n][t];
				cut_fraction_unsel = mt.thinstrength_unsel[n][t];
				patch.man_strength = cut_fraction;
				distribute_cutting(patch, mt.thinselectdiam[n][t], mt.thinselectage[n][t], mt.thinselectpft[n][t], cut_fraction_unsel);
			}
		}
	}
}

/// Harvest of tree individuals by an amount man_strength by calling harvest_wood()
/**	Various carbon accounting for printing purposes is done here.
 *  If clearcut is selected (depending on result from manage_forest()), the individual is killed.
 *  INPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following public members:
 *   - man_strength 				management strength (cutting intensity)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *									Patch (accessed from indiv) public members:
 *   - man_strength					management strength (cutting intensity)
 *   - distributed_cutting			whether patch man_strength has already been distributed to individuals this year
 *									ManagementType (accessed from indiv) public members:
 *   - killgrass_at_cc				whether grass should be killed during clear-cut
 *   - harv_eff_cc					harvest efficiency (fraction removed) for all pfts during clear-cut
 *   - res_outtake_twig_cc			fraction of twigs and branches removed for all pfts during clear-cut
 *   - res_outtake_coarse_root_cc	fraction of coarse roots removed for all pfts during clear-cut
 *   - harv_eff_thin				harvest efficiency (fraction removed) for all pfts during thinning
 *   - res_outtake_twig_thin		fraction of twigs and branches removed for all pfts during thinning
 *   - res_outtake_coarse_root_thin	fraction of coarse roots removed for all pfts during thinning
 *  \param alive					whether individual has survived the first year (needed when function called last day of the year)
 *  \param anpp						individual npp this year (needed when function called last day of the year)
 *  \param pft						reference to a Pft object containing the following public members:
 *   - harv_eff						harvest efficiency (fraction removed) for pft
 *   - res_outtake					fraction of wood residues removed for pft during harvest
 *  OUTPUT PARAMETERS
 *  \param killed					whether individual was killed in this function (needed when function called last day of the year)
 *									Patch (accessed from indiv) public members:
 *   - managed_this_year			whether patch has been cut this year
 *   - has_been_cut					whether patch has ever been cut
 *									Patchpft (accessed from indiv) public members:
 *   - cmass_wood_clearcut			removed wood C biomass of trees harvested during clear-cut (kgC/m2)
 *   - cmass_harv_killed			total C biomass of trees harvested during clear-cut (kgC/m2)
 */
void harvest_forest(Individual& indiv, Pft& pft, bool alive, double anpp, bool& killed) {

	Patch& patch = indiv.vegetation.patch;
	Patchpft& ppft = patch.pft[indiv.pft.id];

	double man_strength = patch.man_strength;
	if(patch.distributed_cutting)
		man_strength = indiv.man_strength;
	ManagementType& mt = indiv.vegetation.patch.stand.get_current_management();

	// Whether to kill grass at clear-cut.
	bool killgrass = mt.killgrass_at_cc && patch.man_strength == 1.00;

	if (man_strength > 0.00 && pft.lifeform == TREE || killgrass && pft.lifeform == GRASS) {

		// Default forestry harvest parameters: pft values
		double harv_eff_wood_harvest = pft.harv_eff;				// Default pft value for trees in global/europe.ins: 0.9
		double res_outtake_twig_wood_harvest = pft.res_outtake;		// Default pft value for trees in global/europe.ins: 0.4
		double res_outtake_coarse_root_wood_harvest = 0.1;			// Using default value of 0.1,

		// Use ManagementType values, if defined
		if (patch.man_strength == 1.00) {
			// the whole cohort is cut (clear-cut)
			if(mt.harv_eff_cc != -1.0)
				harv_eff_wood_harvest = mt.harv_eff_cc;
			if(mt.res_outtake_twig_cc != -1.0)
				res_outtake_twig_wood_harvest = mt.res_outtake_twig_cc;
			if(mt.res_outtake_coarse_root_cc != -1.0)
				res_outtake_coarse_root_wood_harvest = mt.res_outtake_coarse_root_cc;
			ppft.cmass_wood_clearcut += man_strength * check_harvest_cmass(indiv, true);
		}
		else {
			// part of the cohort is cut (thinning)
			if(mt.harv_eff_thin != -1.0)
				harv_eff_wood_harvest = mt.harv_eff_thin;
			if(mt.res_outtake_twig_thin != -1.0)
				res_outtake_twig_wood_harvest = mt.res_outtake_twig_thin;
			if(mt.res_outtake_coarse_root_thin != -1.0)
				res_outtake_coarse_root_wood_harvest = mt.res_outtake_coarse_root_thin;
		}

		// Non-commercial harvest, set to true during first thinning in a clearcut management scheme
		if(patch.harvest_to_litter) {
			harv_eff_wood_harvest = 0.0;
			res_outtake_twig_wood_harvest = 0.0;
			res_outtake_coarse_root_wood_harvest = 0.0;
		}
		ppft.cmass_harv_killed += man_strength * indiv.ccont();

		harvest_wood(indiv, man_strength, harv_eff_wood_harvest, res_outtake_twig_wood_harvest, res_outtake_coarse_root_wood_harvest);
		// Kill grass at clear-cut and send biomass to litter if mt.killgrass_at_cc == true
		if(pft.lifeform == GRASS)
			indiv.kill();

		indiv.densindiv *= (1.0 - man_strength);

		if (negligible(indiv.densindiv)) {
			// C balance accounting when anpp is non-zero. This only happens when this function is called on another day
			// than the first day of the year.
			if(indiv.alive && pft.lifeform == TREE) {
				if(anpp > 0.0) {
					ppft.cmass_litter_sap += anpp;
				}
				else {
					patch.fluxes.report_flux(Fluxes::HARVESTC, anpp);
				}
			}
			indiv.vegetation.killobj();
			killed = true;
		}
		else {
			allometry(indiv);
		}

		patch.managed_this_year = true;		
		patch.has_been_cut = true;
	}
}

/// Harvest function for pasture, representing grazing (previous year).
/**  Function for balancing carbon and nitrogen fluxes from last year's growth
 *  A fraction of leaves is harvested (pft.harv_eff) and returned as acflux_harvest
 *  This represents grazing minus return as manure.
 *  The rest is handled like natural grass in turnover().
 *  Called from growth() last day of the year for normal harvest/grazing.
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  INPUT PARAMETERS
 *  \param pft						reference to a Pft containing the following public members:
 *   - harv_eff    					harvest efficiency (fraction removed) for pft
 *   - harvest_slow_frac			fraction of harvested products that goes to long-lived products
 *   - res_outtake    				fraction of residue outtake at harvest
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members copied to and from the corresponding
 *									variables of an Individual:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *									,and the following public members copied to and from the corresponding variables of
 *									a Patchpft containing an Individual:
 *   - harvested_products_slow		harvest product pool (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen product pool (kgN/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members added to the corresponding variables
 *									of a Patch containing an individual
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgN/m2)
 */
void harvest_pasture(Harvest_CN& i, Pft& pft, bool alive) {

	double harvest;

	// harvest of leaves (grazing)

	// Carbon:
	harvest = pft.harv_eff * i.cmass_leaf;

	if (ifslowharvestpool) {
		i.cmass_harvested_products_slow += harvest * pft.harvest_slow_frac;
		harvest = harvest * (1 - pft.harvest_slow_frac);
	}
	if (alive)
		i.acflux_harvest += harvest;
	i.cmass_leaf -= harvest;

	// Nitrogen
	harvest = pft.harv_eff * i.nmass_leaf;

	if (ifslowharvestpool) {
		i.nmass_harvested_products_slow += harvest * pft.harvest_slow_frac;
		harvest = harvest * (1 - pft.harvest_slow_frac);
	}
	i.nmass_leaf -= harvest;

	// Reduced removal of N relative to C during grazing.
	double N_harvest_scale = 0.25; // Value that works. Needs to be verified in literature.
	i.anflux_harvest += harvest * N_harvest_scale;
	i.nmass_litter_leaf += harvest * (1.0 - N_harvest_scale);

	if (grassforcrop && alive) {
		// Carbon:
		double residue_outtake = pft.res_outtake * i.cmass_leaf;	// res_outtake currently set to 0.0,
		i.acflux_harvest += residue_outtake;				// could be used for burning
		i.cmass_leaf -= residue_outtake;

		// Nitrogen:
		residue_outtake = pft.res_outtake * i.nmass_leaf;
		i.anflux_harvest += residue_outtake;
		i.nmass_leaf -= residue_outtake;
	}
}

/// Harvest function for pasture, representing grazing (previous year).
/*  Function for balancing carbon and nitrogen fluxes from last year's growth
 *  A fraction of leaves is harvested (pft.harv_eff) and returned as acflux_harvest
 *  This represents grazing minus return as manure.
 *  The rest is handled like natural grass in turnover().
 *  Called from growth() last day of the year for normal harvest/grazing.
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_pasture function.
 *  After the execution of the main harvest_pasture function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *									Patchpft (accessed from indiv) public members:
 *   - harvested_products_slow		harvest product pool (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen product pool (kgN/m2)
 *									Patch (accessed from indiv) public members:
 *   - Fluxes::HARVESTC				harvest flux to atmosphere (kgC/m2)
 *   - Fluxes::HARVESTN   			harvest nitrogen flux out of system (kgN/m2)
 *   								Landcover (accessed from indiv) public members:
 *   - acflux_landuse_change		gridcell-level C flux from harvest associated with landcover change (kgC/m2)
 *   - acflux_landuse_change[lc]	landcover-level C flux from harvest associated with landcover change (kgN/m2)
 */
void harvest_pasture(Individual& indiv, Pft& pft, bool alive, bool lc_change) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	harvest_pasture(indiv_cp, pft, alive);

	indiv_cp.copy_to_indiv(indiv, false, lc_change);

	if(lc_change) {

		Stand& stand = indiv.vegetation.patch.stand;

		stand.get_gridcell().landcover.acflux_landuse_change 
			+= stand.get_gridcell_fraction() * indiv_cp.acflux_harvest /(double)stand.nobj;
		stand.get_gridcell().landcover.anflux_landuse_change 
			+=stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;

		if(stand.lc_origin < NLANDCOVERTYPES) {
			stand.get_gridcell().landcover.acflux_landuse_change_lc[stand.lc_origin] 
				+= stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
			stand.get_gridcell().landcover.anflux_landuse_change_lc[stand.lc_origin] 
				+= stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
		}
	}
}

/// Harvest function for cropland, including true crops, intercrop grass
/**   and pasture grass grown in cropland.
 *  Function for balancing carbon and nitrogen fluxes from this year's harvested carbon and nitrogen.
 *  A fraction of harvestable organs (grass:leaves) is harvested (pft.harv_eff) and returned as acflux_harvest.
 *  A fraction of leaves is removed (pft.res_outtake) and returned as acflux_harvest
 *  The rest, including roots, is returned as litter, leaving NO carbon or nitrogen in living tissue.
 *  Called from growth() last day of the year for old-style harvest/grazing or, alternatively, from crop_growth_daily()
 *	at harvest day (hdate) or last intercrop day (eicdate).
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  This function takes a Harvest_CN struct as an input parameter, copied from an individual and it's associated patchpft and patch.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param isintercropgrass			whether individual is cover crop grass
 *  \param pft						reference to a Pft containing the following public members:
 *   - harv_eff    					harvest efficiency (fraction removed) for pft
 *   - phenology   					leaf phenology (cropgreen, any)
 *   - harvest_slow_frac			fraction of harvested products that goes to long-lived products
 *   - res_outtake    				fraction of residue outtake at harvest
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members copied to and from the corresponding
 *									variables of an Individual:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - cmass_dead_leaf				dead leaf C biomass (kgC/m2)
 *   - cmass_stem					stem pool C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_ho						harvestable organ nitrogen biomass (kgN/m2)
 *   - nmass_agpool					above-ground pool nitrogen biomass (kgN/m2)
 *   - nmass_dead_leaf				dead leaf N biomass (kgN/m2)
 *   - nstore_labile    			labile nitrogen storage (kgN/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgN/m2)
 *									,and the following public members copied to and from the corresponding variables of
 *									a Patchpft containing an Individual:
 *   - litter_leaf    				leaf C litter (kgC/m2)
 *   - litter_root 					root C litter (kgC/m2)
 *   - nmass_litter_leaf 			leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			root nitrogen litter (kgN/m2)
 *   - harvested_products_slow		harvest product pool (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen product pool (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members added to the corresponding variables
 *									of a Patch containing an individual
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgN/m2)
 */
void harvest_crop(Harvest_CN& i, Pft& pft, bool alive, bool isintercropgrass) {

	double residue_outtake, harvest;

	if (pft.phenology==CROPGREEN) {

	// all root carbon and nitrogen goes to litter
		if (i.cmass_root > 0.0)
			i.cmass_litter_root += i.cmass_root;
		i.cmass_root = 0.0;

		if (i.nmass_root > 0.0)
			i.nmass_litter_root += i.nmass_root;
		if (i.nstore_labile > 0.0)
			i.nmass_litter_root += i.nstore_labile;
		if (i.nstore_longterm > 0.0)
			i.nmass_litter_root += i.nstore_longterm;
		i.nmass_root = 0.0;
		i.nstore_labile = 0.0;
		i.nstore_longterm = 0.0;

		// harvest of harvestable organs
		// Carbon:
		if (i.cmass_ho > 0.0) {
			// harvested products
			harvest = pft.harv_eff * i.cmass_ho;

			// not removed harvestable organs are put into litter
			if (pft.aboveground_ho)
				i.cmass_litter_leaf += (i.cmass_ho - harvest);
			else
				i.cmass_litter_root += (i.cmass_ho - harvest);

			// harvested products not consumed (oxidised) this year put into harvested_products_slow
			if (ifslowharvestpool) {
				i.cmass_harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}

			// harvested products consumed (oxidised) this year put into acflux_harvest
			i.acflux_harvest += harvest;
		}
		i.cmass_ho = 0.0;

		// Nitrogen:
		if (i.nmass_ho > 0.0) {

			// harvested products
			harvest = pft.harv_eff * i.nmass_ho;

			// not removed harvestable organs are put into litter
			if (pft.aboveground_ho)
				i.nmass_litter_leaf += (i.nmass_ho - harvest);
			else
				i.nmass_litter_root += (i.nmass_ho - harvest);

			// harvested products not consumed this year put into harvested_products_slow_nmass
			if (ifslowharvestpool) {
				i.nmass_harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}

			// harvested products consumed this year put into anflux_harvest
			i.anflux_harvest += harvest;
		}
		i.nmass_ho = 0.0;

		// residues

		// Carbon
		if ((i.cmass_leaf + i.cmass_agpool + i.cmass_dead_leaf + i.cmass_stem) > 0.0) {

			// removed residues are oxidised
			residue_outtake = pft.res_outtake * (i.cmass_leaf + i.cmass_agpool + i.cmass_dead_leaf + i.cmass_stem);
			i.acflux_harvest += residue_outtake;

			// not removed residues are put into litter
			i.cmass_litter_leaf += i.cmass_leaf + i.cmass_agpool + i.cmass_dead_leaf + i.cmass_stem - residue_outtake;
		}
		i.cmass_leaf = 0.0;
		i.cmass_agpool = 0.0;
		i.cmass_dead_leaf = 0.0;
		i.cmass_stem = 0.0;

		// Nitrogen:
		if ((i.nmass_leaf + i.nmass_agpool + i.nmass_dead_leaf) > 0.0) {

			// removed residues are oxidised
			residue_outtake = pft.res_outtake * (i.nmass_leaf + i.nmass_agpool + i.nmass_dead_leaf);
			i.nmass_litter_leaf += i.nmass_leaf + i.nmass_agpool + i.nmass_dead_leaf - residue_outtake;

			// not removed residues are put into litter
			i.anflux_harvest += residue_outtake;
		}
		i.nmass_leaf = 0.0;
		i.nmass_agpool = 0.0;
		i.nmass_dead_leaf = 0.0;
	}
	else if (pft.phenology == ANY) {

		// Intercrop grass
		if (isintercropgrass) {

			// roots

			// all root carbon and nitrogen goes to litter
			if (i.cmass_root > 0.0)
				i.cmass_litter_root += i.cmass_root;
			if (i.nmass_root > 0.0)
				i.nmass_litter_root += i.nmass_root;
			if (i.nstore_labile > 0.0)
				i.nmass_litter_root += i.nstore_labile;
			if (i.nstore_longterm > 0.0)
				i.nmass_litter_root += i.nstore_longterm;

			i.cmass_root = 0.0;
			i.nmass_root = 0.0;
			i.nstore_labile = 0.0;
			i.nstore_longterm = 0.0;


			// leaves

			// Carbon:
			if (i.cmass_leaf > 0.0) {

				// Harvest/Grazing of leaves:
				harvest = pft.harv_eff_ic * i.cmass_leaf;	// currently no harvest of intercrtop grass

				// not removed grass is put into litter
				i.cmass_litter_leaf += i.cmass_leaf - harvest;

				if (ifslowharvestpool) {
					i.cmass_harvested_products_slow += harvest * pft.harvest_slow_frac; // no slow harvest for grass
					harvest = harvest * (1 - pft.harvest_slow_frac);
				}

				i.acflux_harvest += harvest;
			}
			i.cmass_leaf = 0.0;
			i.cmass_ho = 0.0;
			i.cmass_agpool = 0.0;

			// Nitrogen:
			if (i.nmass_leaf > 0.0) {

				// Harvest/Grazing of leaves:
				harvest = pft.harv_eff_ic * i.nmass_leaf;	// currently no harvest of intercrtop grass

				// not removed grass is put into litter
				i.nmass_litter_leaf += i.nmass_leaf - harvest;

				if (ifslowharvestpool) {
					i.nmass_harvested_products_slow += harvest * pft.harvest_slow_frac; // no slow harvest for grass
					harvest = harvest * (1 - pft.harvest_slow_frac);
				}

				i.anflux_harvest += harvest;
			}
			i.nmass_leaf = 0.0;
			i.nmass_ho = 0.0;
			i.nmass_agpool = 0.0;

		}
		else {	// pasture grass

			// harvest of leaves (grazing)

			// Carbon:
			harvest = pft.harv_eff * i.cmass_leaf;

			if (ifslowharvestpool) {
				i.cmass_harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}
			if (alive)
				i.acflux_harvest += harvest;
			i.cmass_leaf -= harvest;

			i.cmass_ho = 0.0;
			i.cmass_agpool = 0.0;

			// Nitrogen:
			// Reduced removal of N relative to C during grazing.
			double N_harvest_scale = 0.25; // Value that works. Needs to be verified in literature.
			harvest = pft.harv_eff * i.nmass_leaf * N_harvest_scale;

			if (ifslowharvestpool) {
				i.nmass_harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}
			i.anflux_harvest += harvest;
			i.nmass_leaf -= harvest;

			i.nmass_ho=0.0;
			i.nmass_agpool=0.0;
		}
	}
}

/// Harvest function for cropland, including true crops, intercrop grass
/**   and pasture grass grown in cropland.
 *  Function for balancing carbon and nitrogen fluxes from this year's harvested carbon and nitrogen.
 *  A fraction of harvestable organs (grass:leaves) is harvested (pft.harv_eff) and returned as acflux_harvest.
 *  A fraction of leaves is removed (pft.res_outtake) and returned as acflux_harvest
 *  The rest, including roots, is returned as litter, leaving NO carbon or nitrogen in living tissue.
 *  Called from growth() last day of the year for old-style harvest/grazing or, alternatively, from crop_growth_daily()
 *	at harvest day (hdate) or last intercrop day (eicdate).
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_crop() function.
 *  After the execution of the main harvest_crop function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param isintercropgrass			whether individual is cover crop grass
 *  \param harvest_grsC				whether harvest daily carbon values are harvested
 *
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - cmass_dead_leaf				dead leaf C biomass (kgC/m2)
 *   - cmass_stem					tem pool C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_ho						harvestable organ nitrogen biomass (kgN/m2)
 *   - nmass_agpool					above-ground pool nitrogen biomass (kgN/m2)
 *   - nmass_dead_leaf				dead leaf N biomass (kgN/m2)
 *   - nstore_labile    			labile nitrogen storage (kgN/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgN/m2)
 *									Patchpft (accessed from indiv) public members:
 *   - litter_leaf    				leaf C litter (kgC/m2)
 *   - litter_root 					root C litter (kgC/m2)
 *   - nmass_litter_leaf 			leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			root nitrogen litter (kgN/m2)
 *   - harvested_products_slow		harvest product pool (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen product pool (kgN/m2)
 *									Patch (accessed from indiv) public members:
 *   - Fluxes::HARVESTC				harvest flux to atmosphere (kgC/m2)
 *   - Fluxes::HARVESTN   			harvest nitrogen flux out of system (kgC/m2)
 *   								Landcover (accessed from indiv) public members:
 *   - acflux_landuse_change		gridcell-level C flux from harvest associated with landcover change (kgC/m2)
 *   - acflux_landuse_change[lc]	landcover-level C flux from harvest associated with landcover change (kgC/m2)
 */
void harvest_crop(Individual& indiv, Pft& pft, bool alive, bool isintercropgrass, bool harvest_grsC) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv, harvest_grsC);

	harvest_crop(indiv_cp, pft, alive, isintercropgrass);

	indiv_cp.copy_to_indiv(indiv, harvest_grsC);

}

/// Transfers all carbon and nitrogen from living tissue to litter
/** Mainly used at land cover change when remaining vegetation after harvest (grass) is
 *   killed by tillage, following an optional burning. No C goes to product pool.
 *
 *  This function takes a Harvest_CN struct as an input parameter, copied from an individual and it's associated patchpft and patch.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param isintercropgrass			whether individual is cover crop grass
 *  \param burn						whether above-ground vegetation C & N is sent to the atmosphere rather than to litter	
 *  \param pft						reference to a Pft containing the following public members:
 *   - landcover    				specifies type of landcover pft is allowed to grow in
 *   - aboveground_ho   			whether harvestable organs are above ground
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members copied to and from the corresponding
 *									variables of an Individual:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgNm2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgN/m2)
 *   - nmass_ho						harvestable organ nitrogen biomass (kgN/m2)
 *   - nmass_agpool					above-ground pool nitrogen biomass (kgN/m2)
 *   - nstore_labile    			labile nitrogen storage (kgN/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgN/m2)
 *									,and the following public members copied to and from the corresponding variables of
 *									a Patchpft containing an Individual:
 *   - litter_leaf    				leaf C litter (kgC/m2)
 *   - litter_root 					root C litter (kgC/m2)
 *   - nmass_litter_leaf 			leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			root nitrogen litter (kgN/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following public members added to the corresponding variables
 *									of a Patch containing an individual
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgN/m2)
 */
void kill_remaining_vegetation(Harvest_CN& cp, Pft& pft, bool alive, bool istruecrop_or_intercropgrass, bool burn) {


	if (alive || istruecrop_or_intercropgrass)  {
		cp.cmass_litter_root += cp.cmass_root;
		cp.acflux_harvest_tolitter += cp.cmass_root;

		if (burn) {
			cp.acflux_harvest += cp.cmass_leaf;
			cp.acflux_harvest += cp.cmass_sap;
			cp.acflux_harvest += cp.cmass_heart - cp.cmass_debt;
		}
		else {
			cp.cmass_litter_leaf += cp.cmass_leaf;
			cp.cmass_litter_sap += cp.cmass_sap;
			cp.cmass_litter_heart += cp.cmass_heart - cp.cmass_debt;
			cp.acflux_harvest_tolitter += cp.cmass_leaf + cp.cmass_sap + cp.cmass_heart - cp.cmass_debt;
		}
	}

	cp.nmass_litter_root += cp.nmass_root;
	cp.nmass_litter_root += cp.nstore_longterm;
	cp.nmass_litter_root += cp.nstore_labile;

	if (burn) {
		cp.anflux_harvest += cp.nmass_leaf;
		cp.anflux_harvest += cp.nmass_sap;
		cp.anflux_harvest += cp.nmass_heart;
	}
	else {
		cp.nmass_litter_leaf += cp.nmass_leaf;
		cp.nmass_litter_sap += cp.nmass_sap;
		cp.nmass_litter_heart += cp.nmass_heart;
	}

	if (pft.landcover == CROPLAND) {
		if (pft.aboveground_ho) {
			if (burn) {
				cp.acflux_harvest += cp.cmass_ho;
				cp.anflux_harvest += cp.nmass_ho;
			}
			else {
				cp.cmass_litter_leaf += cp.cmass_ho;
				cp.nmass_litter_leaf += cp.nmass_ho;
			}
		}
		else {
			cp.cmass_litter_root += cp.cmass_ho;
			cp.nmass_litter_root += cp.nmass_ho;
		}

		if (burn) {
			cp.acflux_harvest += cp.cmass_agpool;
			cp.anflux_harvest += cp.nmass_agpool;
		}
		else {
			cp.cmass_litter_leaf += cp.cmass_agpool;
			cp.nmass_litter_leaf += cp.nmass_agpool;
		}
	}

	cp.cmass_leaf = 0.0;
	cp.cmass_root = 0.0;
	cp.cmass_sap = 0.0;
	cp.cmass_heart = 0.0;
	cp.cmass_debt = 0.0;
	cp.cmass_ho = 0.0;
	cp.cmass_agpool = 0.0;
	cp.nmass_leaf = 0.0;
	cp.nmass_root = 0.0;
	cp.nstore_longterm = 0.0;
	cp.nstore_labile = 0.0;
	cp.nmass_sap = 0.0;
	cp.nmass_heart = 0.0;
	cp.nmass_ho = 0.0;
	cp.nmass_agpool = 0.0;

}

/// Transfers all carbon and nitrogen from living tissue to litter
/** Mainly used at land cover change when remaining vegetation after harvest (grass) is
 *   killed by tillage, following an optional burning.
 *
 *  This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_crop() function.
 *  After the execution of the main harvest_crop function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param lc_change				whether to save harvest in gridcell-level lc struct
 *  \param burn						whether above-ground vegetation C & N is sent to the atmosphere
 *								     rather than to litter
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgNm2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgN/m2)
 *   - nmass_ho						harvestable organ nitrogen biomass (kgN/m2)
 *   - nmass_agpool					above-ground pool nitrogen biomass (kgN/m2)
 *   - nstore_labile    			labile nitrogen storage (kgN/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgN/m2)
 *									Patchpft (accessed from indiv) public members:
 *   - litter_leaf    				leaf C litter (kgC/m2)
 *   - litter_root 					root C litter (kgC/m2)
 *   - nmass_litter_leaf 			leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			root nitrogen litter (kgN/m2)
 *									Patch (accessed from indiv) public members:
 *   - Fluxes::HARVESTC				harvest flux to atmosphere (kgC/m2)
 *   - Fluxes::HARVESTN   			harvest nitrogen flux out of system (kgC/m2)
 *   								Landcover (accessed from indiv) public members:
 *   - acflux_landuse_change		gridcell-level C flux from harvest associated with landcover change (kgC/m2)
 *   - acflux_landuse_change[lc]	landcover-level C flux from harvest associated with landcover change (kgC/m2)
 */
void kill_remaining_vegetation(Individual& indiv, bool burn, bool lc_change) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	kill_remaining_vegetation(indiv_cp, indiv.pft, indiv.alive, indiv.istruecrop_or_intercropgrass(), burn);

	indiv_cp.copy_to_indiv(indiv, false, lc_change);

	if (burn && lc_change) {
		Stand& stand = indiv.vegetation.patch.stand;
		Landcover& lc = stand.get_gridcell().landcover;
		lc.acflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
		lc.anflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
		if(stand.lc_origin < NLANDCOVERTYPES) {
			lc.acflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
			lc.anflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
		}
	}

}

/// Scaling of last year's or harvest day individual carbon and nitrogen member values in stands that have increased their area fraction this year.
/** Called immediately before harvest functions in growth() or allocation_crop_daily().
 *  INPUT PARAMETERS
 *  \param scale_grsC				whether to scale growing-season biomass variables
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - cmass_dead_leaf				dead leaf C biomass (kgC/m2)
 *   - cmass_stem					stem pool C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - grs_cmass_leaf 				growing-season leaf C biomass (kgC/m2)
 *   - grs_cmass_root				growing-season fine root C biomass (kgC/m2)
 *   - grs_cmass_ho					growing-season harvestable organ C biomass (kgC/m2)
 *   - grs_cmass_agpool				growing-season above-ground pool C biomass (kgC/m2)
 *   - grs_cmass_dead_leaf			growing-season dead leaf C biomass (kgC/m2)
 *   - grs_cmass_stem				growing-season stem pool C biomass (kgC/m2)
 *   - grs_cmass_leaf_luc 			growing-season leaf C biomass on day of land-use change (kgC/m2)
 *   - grs_cmass_root_luc			growing-season fine root C biomass on day of land-use change (kgC/m2)
 *   - grs_cmass_ho_luc				growing-season harvestable organ C biomass on day of land-use change (kgC/m2)
 *   - grs_cmass_agpool_luc			growing-season above-ground pool C biomass on day of land-use change (kgC/m2)
 *   - grs_cmass_dead_leaf_luc		growing-season dead leaf C biomass on day of land-use change (kgC/m2)
 *   - grs_cmass_stem_luc			growing-season stem pool C biomass on day of land-use change (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - param nmass_ho				harvestable organ nitrogen biomass (kgN/m2)
 *   - param nmass_agpool			above-ground pool nitrogen biomass (kgN/m2)
 *   - nmass_dead_leaf				dead leaf N biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgNm2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgN/m2)
 *   - nstore_labile    			labile nitrogen storage (kgN/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgN/m2)
 *   - nmass_leaf_luc 				leaf nitrogen biomass on day of land-use change (kgN/m2)
 *   - nmass_root_luc 				fine root nitrogen biomass on day of land-use change (kgN/m2)
 *   - nmass_ho_luc					harvestable organ nitrogen biomass on day of land-use change (kgN/m2)
 *   - nmass_agpool_luc				above-ground pool nitrogen biomass on day of land-use change (kgN/m2)
 *   - nmass_dead_leaf_luc			dead leaf N biomass on day of land-use change (kgN/m2)
 *   - nmass_sap_luc   				sapwood nitrogen biomass on day of land-use change (kgNm2)
 *   - nmass_heart_luc    			heartwood nitrogen biomass on day of land-use change (kgN/m2)
 *   - nstore_labile_luc    		labile nitrogen storage on day of land-use change (kgN/m2)
 *   - nstore_longterm_luc    		longterm nitrogen storage on day of land-use change (kgN/m2)
 *									Stand (accessed from indiv) public members:
 *   - scale_LC_change		   		scaling factor for stands that have grown in area this year (old fraction/new fraction)
 */
void scale_indiv(Individual& indiv, bool scale_grsC) {

	Stand& stand = indiv.vegetation.patch.stand;
	Gridcell& gridcell = stand.get_gridcell();

	if (stand.scale_LC_change >= 1.0) {
		return;
	}

	// Scale individual's C and N mass in stands that have increased in area
	// this year by (old area/new area):
	double scale = stand.scale_LC_change;

	if (scale_grsC) {

		if (indiv.pft.landcover == CROPLAND) {

			if (indiv.has_daily_turnover()) {

				indiv.cropindiv->grs_cmass_leaf -= indiv.cropindiv->grs_cmass_leaf_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_root -= indiv.cropindiv->grs_cmass_root_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_ho -= indiv.cropindiv->grs_cmass_ho_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_agpool -= indiv.cropindiv->grs_cmass_agpool_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_dead_leaf -= indiv.cropindiv->grs_cmass_dead_leaf_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_stem -= indiv.cropindiv->grs_cmass_stem_luc * (1.0 - scale);

				indiv.check_C_mass();
			}
			else {
				indiv.cropindiv->grs_cmass_leaf *= scale;
				indiv.cropindiv->grs_cmass_root *= scale;
				indiv.cropindiv->grs_cmass_ho *= scale;
				indiv.cropindiv->grs_cmass_agpool *= scale;
				indiv.cropindiv->grs_cmass_plant *= scale;	//grs_cmass_plant not used
				indiv.cropindiv->grs_cmass_dead_leaf *= scale;
				indiv.cropindiv->grs_cmass_stem *= scale;
			}
		}
	}
	else {

		indiv.cmass_root *= scale;
		indiv.cmass_leaf *= scale;
		indiv.cmass_heart *= scale;
		indiv.cmass_sap *= scale;
		indiv.cmass_debt *= scale;

		if (indiv.pft.landcover == CROPLAND) {
			indiv.cropindiv->cmass_agpool *= scale;
			indiv.cropindiv->cmass_ho *= scale;
		}
	}

	// Deduct individual N present day 0 this year in stands that have increased in area this year,
	// scaled by (1 - old area/new area):
	indiv.nmass_root = indiv.nmass_root - indiv.nmass_root_luc * (1.0 - scale);
	indiv.nmass_leaf = indiv.nmass_leaf - indiv.nmass_leaf_luc * (1.0 - scale);
	indiv.nmass_heart = indiv.nmass_heart - indiv.nmass_heart_luc * (1.0 - scale);
	indiv.nmass_sap = indiv.nmass_sap - indiv.nmass_sap_luc * (1.0 - scale);

	if (indiv.pft.landcover == CROPLAND) {
		indiv.cropindiv->nmass_agpool = indiv.cropindiv->nmass_agpool - indiv.cropindiv->nmass_agpool_luc * (1.0 - scale);
		indiv.cropindiv->nmass_ho = indiv.cropindiv->nmass_ho - indiv.cropindiv->nmass_ho_luc * (1.0 - scale);
		indiv.cropindiv->nmass_dead_leaf =indiv.cropindiv->nmass_dead_leaf - indiv.cropindiv->nmass_dead_leaf_luc * (1.0 - scale);
	}

	if (indiv.nstore_labile > indiv.nstore_labile_luc * (1.0 - scale)) {
		indiv.nstore_labile -= indiv.nstore_labile_luc * (1.0 - scale);
	}
	else {
		indiv.nstore_longterm -= indiv.nstore_labile_luc * (1.0 - scale);
	}
	indiv.nstore_longterm = indiv.nstore_longterm - indiv.nstore_longterm_luc * (1.0 - scale);

	indiv.check_N_mass();
}

/// Yearly function for harvest of cropland and pasture stands that have yearly allocation, turnover and gridcell.expand_to_new_stand[lc] = false.
/** Should only be called from growth().
//  Harvest functions are preceded by rescaling of living C.
//  Only affects natural stands if gridcell.expand_to_new_stand[NATURAL] is false.
 */
bool harvest_year(Individual& indiv) {

	Stand& stand = indiv.vegetation.patch.stand;
	Landcover& landcover = stand.get_gridcell().landcover;
	bool killed = false;

	// Reduce individual's C and N mass in stands that have increased in area this year:
	if (landcover.updated && !indiv.has_daily_turnover()) {
		scale_indiv(indiv, false);
	}

	if (stand.landcover == CROPLAND) {
		if (!indiv.has_daily_turnover())
			harvest_crop(indiv, indiv.pft, indiv.alive, indiv.cropindiv->isintercropgrass, false);
	}
	else if (stand.landcover == PASTURE) {
		harvest_pasture(indiv, indiv.pft, indiv.alive);
	}

	return killed;
}

/// Yield function for true crops and inter-crop grass.
/** INPUT PARAMETERS
 *  \param indiv					reference to an Individual
 *									Pft (accessed from indiv) public members:
 *   - phenology 					leaf phenology (cropgreen, any)
 *   - harv_eff 					harvest efficiency of crop
 *   - harv_eff_ic 					harvest efficiency of inter-crop grass
 *									cropindiv_struct (accessed from indiv) public members:struct containing the following
 *									public members:
 *   - ycmass_leaf 					daily updated leaf C biomass (kgC/m2), reset at day 0
 *   - harv_cmass_leaf 				year's harvested leaf C biomass (kgC/m2)
 *   - ycmass_ho 					daily updated harvestable organ C biomass (kgC/m2), reset at day 0
 *   - harv_cmass_ho 				year's harvested harvestable organ C biomass (kgC/m2)
 *   - cmass_ho_harvest[2] 			harvestable organ C biomass at the last two harvest events this year (kgC/m2)	
 *  OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual
 *									cropindiv_struct (accessed from indiv) public members:struct containing the following
 *									public members:
 *   - yield 						dry weight crop yield grown this year (kgC/m2)
 *   - harv_yield 					dry weight crop yield harvested this year (kgC/m2)
 *   - yield_harvest[2] 			dry weight crop yield at the last two harvest events this year (kgC/m2)
 */
void yield_crop(Individual& indiv) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());

	if (indiv.pft.phenology == ANY) {			// grass intercrop yield

		/* Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in
		 * harvest_crop (different years)
		 */
		if (cropindiv.ycmass_leaf > 0.0) {
			cropindiv.yield = cropindiv.ycmass_leaf * indiv.pft.harv_eff_ic * 2.0;
		}
		else {
			cropindiv.yield = 0.0;
		}

		// Yield dry wieght of actually harvest products this year; NB as above
		if (cropindiv.harv_cmass_leaf > 0.0) {
			cropindiv.harv_yield = cropindiv.harv_cmass_leaf * indiv.pft.harv_eff_ic * 2.0;
		}
		else {
			cropindiv.harv_yield = 0.0;
		}
	}
	else if (indiv.pft.phenology == CROPGREEN) {		//true crop yield

		/* Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in
		 * harvest_crop (different years)
		*/
		if (cropindiv.ycmass_ho > 0.0) {
			cropindiv.yield = cropindiv.ycmass_ho * indiv.pft.harv_eff * 2.0;// Should be /0.446 instead
		}
		else {
			cropindiv.yield = 0.0;
		}

		// Yield dry wieght of actually harvest products this year; NB as above
		if (cropindiv.harv_cmass_ho > 0.0) {
			cropindiv.harv_yield=cropindiv.harv_cmass_ho * indiv.pft.harv_eff * 2.0;
		}
		else {
			cropindiv.harv_yield = 0.0;
		}

		// Yield dry wieght of actually harvest products this year; NB as above
		for (int i=0;i<2;i++) {
			if (cropindiv.cmass_ho_harvest[i] > 0.0) {
				cropindiv.yield_harvest[i] = cropindiv.cmass_ho_harvest[i] * indiv.pft.harv_eff * 2.0;
			}
			else {
				cropindiv.yield_harvest[i]=0.0;
			}
		}
	}
}

/// Yield function for pasture grass grown in cropland landcover
/** INPUT PARAMETERS
 *  \param indiv					reference to an Individual
 *  \param cmass_leaf_inc			increase in leaf C this year (kgC/m2)
 *									Pft (accessed from indiv) public members:	
 *   - harv_eff 					harvest efficiency of grass
 *  OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual
 *									cropindiv_struct (accessed from indiv) public members:struct containing the following
 *									public members:
 *   - yield 						dry weight grass yield grown this year (kgC/m2)
 *   - harv_yield 					dry weight grass yield harvested this year (kgC/m2)
 */
void yield_pasture(Individual& indiv, double cmass_leaf_inc) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());

	// Normal CC3G/CC4G stand growth (Pasture)

	// OK if turnover_leaf==1.0, else (cmass_leaf+cmass_leaf_inc)*indiv.pft.harv_eff*2.0
	cropindiv.yield = max(0.0, cmass_leaf_inc) * indiv.pft.harv_eff * 2.0;
	// Although no specified harvest date, harv_yield is set for compatibility.
	cropindiv.harv_yield = cropindiv.yield;
}

/// Function that determines amount of nitrogen applied today. Crop-specific, Pft-based.
/** INPUT PARAMETERS
 *  \param patch					reference to a Patch
 *									Pft (accessed from looping through pftlist) public members:
 *   - phenology 					leaf phenology (cropgreen, any)
 *   - N_appfert 					nitrogen fertilization this year (kgN/m2)
 *   - fertrate[2] 					how much of the fertiliser that is applied at the two fertilisation events after sowing
 *   - fert_stages[2] 				development stage at the two fertilisation events after sowing
 *									cropphen_struct (accessed from patch and pft) public members:
 *   - growingseason 				whether inside crop/intercrop grass growing period
 *   - fertilised[3] 				whether this field was fertilised or not at the three available fertilisation occasions
 *   - dev_stage					crop development stage (0-2)
 *									Gridcellpft (accessed from looping through pftlist) public members:
 *   - Nfert_read					total nitrogen fertilization this year read from text input file (kgN/m2)
 *   - Nfert_man_read				nitrogen fertilization from manure this year read from text input file (kgN/m2)
 *									Gridcellst (accessed from patch) public members:
 *   - nfert						total StandType-level nitrogen fertilization this year read from instruction file or
 *									text input file (kgN/m2)
 *  INPUT/OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - anfert 						year's nitrogen fertilization (kgN/m2)
 *   - Fluxes::NFERT   				fertilisation N flux to vegetation (kgN/m2)
 *   - Fluxes::MANUREC   			C flux to vegetation associated with manure addition (kgC/m2)
 *   - Fluxes::MANUREN   			N flux to vegetation associated with manure addition (kgN/m2)
 *									sompool[] (accessed from patch) containing the following public members:
 *   - nmass 						nitrogen mass in pool (kgN/m2)
 *   - cmass 						C mass in pool (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - dnfert 						nitrogen fertilization today (kgN/m2)
 */
void nfert_crop(Patch& patch) {

	Gridcell& gridcell = patch.stand.get_gridcell();
	StandType& st = stlist[patch.stand.stid];

	patch.dnfert = 0.0;

	pftlist.firstobj();
	// Loop through PFTs
	while(pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Patchpft& patchpft = patch.pft[pft.id];
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if (patch.stand.pft[pft.id].active && pft.phenology == CROPGREEN) {

			cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
			if(!ppftcrop.growingseason) {
				pftlist.nextobj();
				continue;
			}

			double nfert = pft.N_appfert;
			double mineral = 1.0;
			// Use total and manure fertilisation amount from text input file if present:
			if (gridcellpft.Nfert_read >= 0.0) {
				nfert = gridcellpft.Nfert_read;
				if (gridcellpft.Nfert_man_read > 0.0) {
					mineral = 1.0 - gridcellpft.Nfert_man_read;
				}
			}
			/// StandType-level fertilisation input for crops:
			if(gridcell.st[st.id].nfert >= 0.0) {
				nfert = gridcell.st[st.id].nfert;
			}

			if (!ppftcrop.fertilised[0] && ppftcrop.dev_stage > 0.0) {
				// Fertiliser application directly after sowing:
				patch.dnfert = nfert * mineral * (1.0 - pft.fertrate[0] - pft.fertrate[1]);
				patch.fluxes.report_flux(Fluxes::NFERT,nfert * mineral * (1.0 - pft.fertrate[0] - pft.fertrate[1]));
				ppftcrop.fertilised[0] = true;
				if (mineral < 1.0) {
					patch.soil.sompool[SOILMETA].nmass += nfert * (1.0 - mineral) * 0.5;
					patch.soil.sompool[SOILMETA].cmass += nfert * (1.0 - mineral) * 30.0 * 0.25;
					patch.soil.sompool[SOILSTRUCT].nmass += nfert * (1.0 - mineral) * 0.5;
					patch.soil.sompool[SOILSTRUCT].cmass += nfert * (1.0 - mineral) * 30.0 * 0.75;
					patch.fluxes.report_flux(Fluxes::MANUREC, -nfert * (1.0 - mineral) * 30.0 );
					patch.anfert += nfert * (1.0 - mineral);
					patch.fluxes.report_flux(Fluxes::MANUREN, nfert * (1.0 - mineral));
				}
			}
			else if (!ppftcrop.fertilised[1] && ppftcrop.dev_stage > pft.fert_stages[0]) {
				// Fertiliser application at first fertilisation event after sowing:
				patch.dnfert = nfert * mineral * pft.fertrate[0];
				patch.fluxes.report_flux(Fluxes::NFERT,nfert * mineral * pft.fertrate[0]);
				ppftcrop.fertilised[1] = true;
			}
			else if (!ppftcrop.fertilised[2] && ppftcrop.dev_stage > pft.fert_stages[1]) {
				// Fertiliser application at second fertilisation event after sowing:
				patch.dnfert = nfert * mineral * (pft.fertrate[1]);
				patch.fluxes.report_flux(Fluxes::NFERT,nfert * mineral * pft.fertrate[1]);
				ppftcrop.fertilised[2] = true;
			}
		}
		pftlist.nextobj();
	}
	patch.anfert += patch.dnfert;
}

/// Function that determines amount of phosphorus applied today. Crop-specific, pft-based. General, both StandType-based and Pft-based, calling function nfert_crop().
/** INPUT PARAMETERS
 *  \param patch					reference to a Patch
 *									Stand (accessed from patch) public members:
 *   - landcover		   			type of landcover
 *									Gridcellst (accessed from patch) public members:
 *   - nfert						total StandType-level nitrogen fertilization this year read from instruction file or
 *									text input file (kgN/m2)
 *  INPUT/OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - anfert 						year's nitrogen fertilization (kgN/m2)
 *   - Fluxes::NFERT   				fertilisation N flux to vegetation (kgN/m2)
 *  OUTPUT PARAMETERS
 *  \param patch					reference to a Patch containing the following public members:
 *   - dnfert 						nitrogen fertilization today (kgN/m2)
 */

void pfert_crop(Patch& patch) {

	Gridcell& gridcell = patch.stand.get_gridcell();

	patch.dpfert = 0.0;

	pftlist.firstobj();
	// Loop through PFTs
	while (pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Patchpft& patchpft = patch.pft[pft.id];
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if (patch.stand.pft[pft.id].active && pft.phenology == CROPGREEN) {

			cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
			if (!ppftcrop.growingseason) {
				pftlist.nextobj();
				continue;
			}

			double pfert = pft.P_appfert;
			double mineral = 1.0;
			if (gridcellpft.Pfert_read >= 0.0) {
				pfert = gridcellpft.Pfert_read;
				if (gridcellpft.Pfert_man_read > 0.0) {
					mineral = 1.0 - gridcellpft.Pfert_man_read;
				}

			}
			if (!ppftcrop.fertilised[0] && ppftcrop.dev_stage > 0.0) {
				// Fertiliser application at dev_stage = 0, sowing.
				patch.dpfert = pfert * mineral * (1.0 - pft.fertrate[0] - pft.fertrate[1]);
				patch.fluxes.report_flux(Fluxes::PFERT, pfert * mineral * (1.0 - pft.fertrate[0] - pft.fertrate[1]));
				ppftcrop.fertilised[0] = true;
				if (mineral<1.0) {
					patch.soil.sompool[SOILMETA].pmass += pfert * (1.0 - mineral) * 0.5;
					patch.soil.sompool[SOILMETA].cmass += pfert * (1.0 - mineral) * 610.6792 * 0.25; //Check this, C:P ratio of an SLA of 15.51 m2/kgC (same which leads to a C:N of 30)
					patch.soil.sompool[SOILSTRUCT].pmass += pfert * (1.0 - mineral) * 0.5;
					patch.soil.sompool[SOILSTRUCT].cmass += pfert * (1.0 - mineral) * 610.6792 * 0.75;
					patch.fluxes.report_flux(Fluxes::MANUREC, -pfert * (1.0 - mineral) * 610.6792);
					patch.anfert += pfert * (1.0 - mineral);
					patch.fluxes.report_flux(Fluxes::MANUREP, pfert * (1.0 - mineral));
				}
			}
			else if (!ppftcrop.fertilised[1] && ppftcrop.dev_stage > pft.fert_stages[0]) {
				patch.dpfert = pfert * mineral * pft.fertrate[0];
				patch.fluxes.report_flux(Fluxes::PFERT, pfert * mineral * pft.fertrate[0]);
				ppftcrop.fertilised[1] = true;
			}
			else if (!ppftcrop.fertilised[2] && ppftcrop.dev_stage > pft.fert_stages[1]) {
				patch.dpfert = pfert * mineral * (pft.fertrate[1]);
				patch.fluxes.report_flux(Fluxes::PFERT, pfert * mineral * pft.fertrate[1]);
				ppftcrop.fertilised[2] = true;
			}
		}
		pftlist.nextobj();
	}
	patch.apfert += patch.dpfert;
}

/// Function that determines amount of nitrogen applied today.
void nfert(Patch& patch) {

	Stand& stand = patch.stand;
	StandType& st = stlist[stand.stid];
	Gridcell& gridcell = stand.get_gridcell();

	if(stand.landcover == CROPLAND) {
		nfert_crop(patch);
		return;
	}

	// General code for applying nitrogen to other land cover types, an equal amount every day.
	double nfert;
	if(gridcell.st[st.id].nfert >= 0.0) {	// todo: management type variable (mt.nfert)
		nfert = gridcell.st[st.id].nfert;
	}
	else {
		nfert = 0.0;
	}
	patch.dnfert = nfert / date.year_length();
	patch.anfert += patch.dnfert;
	patch.fluxes.report_flux(Fluxes::NFERT, patch.dnfert);
}


/// Function that determines amount ofphosphorus applied today.
void pfert(Patch& patch) {

	Stand& stand = patch.stand;
	StandType& st = stlist[stand.stid];
	Gridcell& gridcell = stand.get_gridcell();

	if (stand.landcover == CROPLAND) {
		pfert_crop(patch);
		return;
	}

	// General code for applying phosphorus to other land cover types, an equal amount every day.
	double pfert;
	if (gridcell.st[st.id].pfert >= 0.0) {	// todo: management type variable (mt.nfert)
		pfert = gridcell.st[st.id].pfert;
	}
	else {
		pfert = 0.0;
	}
	patch.dpfert = pfert / date.year_length();
	patch.apfert += patch.dpfert;
}


// Decides when to apply forest management rotation status by calling stand.rotate()
/** Sets new forest management variables by calling stand.rotate() on st.mtstartyear[m]
*  \param stand					reference to a Stand containing the following public members:
*   - landcover		   			type of landcover
*   - current_rot		   			current management rotation item
*									Patch (accessed from looping through stand) containing the following public members:
*   - clearcut_this_year			whether patch has been clearcut this year
*									StandType (accessed from Stand) public members:
*   - rotation_wait_for_cc			whether to wait for clearcut before moving to next ManagementType in a forestry rotation
*   - mtstartyear[]				start of the managements in a rotation cycle (calendar year)
*									Rotation (accessed from stand) public members:
*   - nmanagements					number of managements in rotation
*									ManagementType (accessed from Stand) public members:
*   - harvest_system				whether clear-cut or continuous cutting schemes are selected
*  INPUT/OUTPUT PARAMETERS
*  \param stand					reference to a Stand containing the following public members:
*   - nyears_in_rotation		   	number of years passed in current rotation item
*/
void forest_rotation(Stand& stand) {

	StandType& st = stlist[stand.stid];
	ManagementType& mt = stand.get_current_management();

	if (stand.landcover != FOREST || st.rotation.nmanagements < 2)
		return;

	stand.nyears_in_rotation++;

	bool clearcut = false;
	double cmass_wood = 0.0;
	for (unsigned int p = 0; p < stand.nobj; p++) {
		if (stand[p].clearcut_this_year)
			clearcut = true;
		cmass_wood += stand[p].cmass_wood();	// don't wait for clearcut to rotate if no trees are alive
	}

	bool rotate_at_cc = st.rotation_wait_for_cc && mt.harvest_system == "CLEARCUT";

	for (int m = 0; m<st.rotation.nmanagements; m++) {
		bool rotate_now = false;
		if (rotate_at_cc) {
			if (st.mtstartyear[m] > -1 && st.mtstartyear[m] <= date.get_calendar_year() && clearcut
				|| st.mtstartyear[m] == date.get_calendar_year() && !cmass_wood) {
				rotate_now = true;
			}
		}
		else if (st.mtstartyear[m] == date.get_calendar_year()) {
			rotate_now = true;
		}
		if (rotate_now && m != stand.current_rot) {
			stand.rotate(m);
		}
	}
}

// Decides when to apply crop management rotation status by calling stand.rotate()

/// Updates crop rotation status

/** Sets new crop management variables, typically on harvest day
 *  \param stand					reference to a Stand containing the following public members:
 *   - landcover		   			type of landcover
 *   - isrotationday		   		whether crop rotation item is to be updated today
 *   - first_year					simulation year when this stand was created.	   		
 *   - current_rot		   			current management rotation item
 *   - pftid		   				pft id of main crop, updated during rotation	
 *   - multicrop		   			whether double cropping of one crop (e.g. rice) occurs
 *									Rotation (accessed from stand) public members:
 *   - nmanagements					number of managements in rotation
 *   - firstrotyear					first crop rotation year
 *									ManagementType (accessed from Stand) public members:
 *   - fallow						whether grass is grown in fallow
 *  INPUT/OUTPUT PARAMETERS
 *  \param stand					reference to a Stand containing the following public members:
 *   - ndays_in_rotation		   		number of days passed in current rotation item
 *   - infallow		   				whether crop stand is in fallow (with cover crop grass)
 *									Standpft (accessed from stand and pftid) public members:
 *   - sdate_force					sowing date specified in stand type or read from input file
 *   - hdate_force					harvest date specified in stand type or read from input file
 *  OUTPUT PARAMETERS
 *  \param stand					reference to a Stand
 *									cropphen_struct (accessed from looping through patches in the stand) public members:
 *   - bicdate		 				day of beginning of intercropseason
 *   - eicdate		 				day of end of intercropseason
 *   - hdate		 				harvest day
 *   - intercropseason		 		whether inside intercrop crass growing period
 *									Gridcellpft (accessed from stand and pftid) public members:
 *   - sowing_restriction			flag to preclude crop sowing during fallow
 */
void crop_rotation(Stand& stand) {

	if (stand.landcover != CROPLAND) {
		return;
	}

	Rotation& rotation = stlist[stand.stid].rotation;

	stand.ndays_in_rotation++;

	if (rotation.nmanagements < 2 || !stand.isrotationday) {
		return;
	}

	int firstrotyear = rotation.firstrotyear - date.first_calendar_year;
	bool postpone_rotation = false;

	if (date.year < stand.first_year + 3) {
		if ((abs(firstrotyear - date.year) % rotation.nmanagements) != stand.current_rot)
			postpone_rotation = true;
	}

	if (!postpone_rotation) {

		if(stand.infallow) {
			stand.infallow = false;
			stand.get_gridcell().pft[stand.pftid].sowing_restriction = false;
		}

		int old_pftid = stand.pftid;

		stand.rotate();

		for (unsigned int p=0; p<stand.nobj; p++) {

			cropphen_struct& previous = *(stand[p].pft[old_pftid].get_cropphen());
			cropphen_struct& current = *(stand[p].pft[stand.pftid].get_cropphen());

			previous.bicdate = -1;
			if(!previous.intercropseason)
				current.bicdate = stepfromdate(date.day, 15);
			previous.eicdate = -1;
			current.eicdate = -1;
			previous.hdate = -1;
			current.intercropseason = previous.intercropseason;
		}

		// Adds sowing and harvest dates for the second crop in a double cropping system
		if (rotation.multicrop && rotation.nmanagements == 2 && stand.current_rot == 1) {
			if (stand.pft[stand.pftid].sdate_force < 0)
				stand.pft[stand.pftid].sdate_force = stepfromdate(date.day, 10);
			if (stand.pft[stand.pftid].hdate_force < 0) {
				stand.pft[stand.pftid].hdate_force = stepfromdate(stand.pft[old_pftid].sdate_force, -10);
			}
		}

		if(stand.get_current_management().fallow) {
			stand.infallow = true;
			stand.get_gridcell().pft[stand.pftid].sowing_restriction = true;
		}
	}

	stand.isrotationday = false;
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Bellassen, V, Le Maire, G, Dhôte, JF & Viovy, N (2010) Modelling forest management within a global vegetation model
//	 - Part 1: Model Structure and general behaviour. Ecol. Modelling 221: 2458-2474.
// Lagergren, F and Jönsson, A M (2017) Ecosystem model analysis of multi-use forestry in a changing climate, Ecosyst.
//	 Serv. 26: 209-224.
// Lindeskog, M, Smith, B, Lagergren, F, Sycheva, E, Ficko, A, Pretzsch, H, and Rammig, A: Accounting for forest management
//	 in the estimation of forest carbon balance using the dynamic vegetation model LPJ-GUESS (v4.0, r9710): implementation
//   and evaluation of simulations for Europe (2021). Geosci. Model Dev. 14: 6071-6112.
