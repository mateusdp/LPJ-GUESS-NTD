////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file management.h
/// \brief Harvest functions for cropland, managed forest and pasture			
/// \author Mats Lindeskog
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_MANAGEMENT_H
#define LPJ_GUESS_MANAGEMENT_H

// Forward declaration
struct Harvest_CN;

/// Harvest function for cropland, including true crops, intercrop grass
void harvest_crop(Harvest_CN& indiv_cp, Pft& pft, bool alive, bool isintercropgrass);

/// Harvest function for cropland, including true crops, intercrop grass
void harvest_crop(Individual& indiv, Pft& pft, bool alive, bool isintercropgrass, bool harvest_grs);

/// Sets forest management for all stands this year
void manage_forests(Gridcell& gridcell);

/// Sets management strength for individual trees to achieve prescribed tree pft composition
void set_forest_pft_structure(Gridcell& gridcell);

/// Harvest function used for managed forest and for clearing natural vegetation at land use change.
void harvest_wood(Harvest_CN& indiv_cp, double diam,Pft& pft, bool alive, double frac_cut, double harv_eff,
				  double res_outtake_twig = 0.0, double res_outtake_coarse_root = 0.0);

/// Harvest function used for managed forest and for clearing natural vegetation at land use change.
void harvest_wood(Individual& indiv, double frac_cut, double harv_eff, double res_outtake_twig = 0.0,
				  double res_outtake_coarse_root = 0.0, bool lc_change = false);

/// Harvest function for pasture, representing grazing.
void harvest_pasture(Harvest_CN& indiv_cp, Pft& pft, bool alive);

/// Harvest function for pasture, representing grazing.
void harvest_pasture(Individual& indiv, Pft& pft, bool alive, bool lc_change = false);

/// Harvest function for managed forest
void harvest_forest(Individual& indiv, Pft& pft, bool alive, double anpp, bool& killed);

/// Transfers all carbon and nitrogen from living tissue to litter.
void kill_remaining_vegetation(Harvest_CN& indiv_cp, Pft& pft, bool alive, bool istruecrop_or_intercropgrass, bool burn = false);

/// Transfers all carbon and nitrogen from living tissue to litter.
void kill_remaining_vegetation(Individual& indiv, bool burn = false, bool lc_change = false);

/// Scaling of last year's or harvest day individual carbon and nitrogen member values in stands that have increased area fraction this year.
void scale_indiv(Individual& indiv, bool scale_grsC);

/// Yearly function for harvest of all land covers. Should only be called from growth()
bool harvest_year(Individual& indiv);

/// Yield function for true crops and intercrop grass
void yield_crop(Individual& indiv);

/// Yield function for pasture grass grown in cropland landcover
void yield_pasture(Individual& indiv, double cmass_leaf_inc);

/// Determines amount of nitrogen applied today
void nfert(Patch& patch);

/// Updates crop rotation status
void crop_rotation(Stand& stand);

/// Updates forestry rotation status
void forest_rotation(Stand& stand);

/// Sets forest management for patch this year
void manage_forest(Patch& patch);

// Returns harvestable cmass for individual
double check_harvest_cmass(Individual& indiv, bool stem_cmass_only = false, bool to_product_pool = false);

// Returns harvestable cmass for patch
double check_harvest_cmass(Patch& patch, bool stem_cmass_only = false, bool check_selection = false);

// Returns harvestable cmass for stand
double check_harvest_cmass(Stand& stand, bool stem_cmass_only = false, bool check_selection = false);

/// Class for ranking stands by harvestable cmass etc.
class Rank_stands {

	Gridcell& gridcell;
	std::vector<double> cmass_harvest_vect;
	std::vector<double> cmass_harvest_wood_vect;
	std::vector<std::pair<double,int> > sort_vect;

public:

	void sort_data(std::vector<double> v) {

		sort_vect.clear();
		sort_vect.reserve(gridcell.nbr_stands());
		for(unsigned int i=0;i<gridcell.nbr_stands();i++) {
			Stand& stand = gridcell[i];
			sort_vect.push_back(std::make_pair(v[i],i));
		}
		std::sort(sort_vect.begin(), sort_vect.end()); 
	}

	void sort_cmass_harvest() {
		sort_data(cmass_harvest_vect);
	}
	void sort_cmass_harvest_wood() {
		sort_data(cmass_harvest_wood_vect);
	}

	Rank_stands(Gridcell& gridcellX) : gridcell(gridcellX) {

		cmass_harvest_vect.reserve(gridcell.nbr_stands());
		cmass_harvest_wood_vect.reserve(gridcell.nbr_stands());
		for(unsigned int i=0;i<gridcell.nbr_stands();i++) {
			Stand& stand = gridcell[i];
			cmass_harvest_vect.push_back(check_harvest_cmass(stand));
			cmass_harvest_wood_vect.push_back(check_harvest_cmass(stand, true));
		}
		// Make sure sort_vect is filled
		sort_cmass_harvest();
	}

	int get_index(int rank) {
		return sort_vect[rank].second;
	}
	double get_cmass_harvest(int rank) {
		return cmass_harvest_vect[sort_vect[rank].second];
	}
	double get_cmass_harvest_wood(int rank) {
		return cmass_harvest_wood_vect[sort_vect[rank].second];
	}
};

/// Class for ranking patches by harvestable cmass etc.
class Rank_patches {

	Stand& stand;
	std::vector<double> cmass_harvest_vect;
	std::vector<double> cmass_harvest_wood_vect;
	std::vector<std::pair<double,int> > sort_vect;

public:

	void sort_data(std::vector<double> v) {

		sort_vect.clear();
		sort_vect.reserve(stand.nobj);
		for(unsigned int i=0;i<stand.nobj;i++) {
			Patch& patch = stand[i];
			sort_vect.push_back(std::make_pair(v[i],i));
		}
		std::sort(sort_vect.begin(), sort_vect.end()); 
	}

	void sort_cmass_harvest() {
		sort_data(cmass_harvest_vect);
	}
	void sort_cmass_harvest_wood() {
		sort_data(cmass_harvest_wood_vect);
	}

	Rank_patches(Stand& standX) : stand(standX) {

		cmass_harvest_vect.reserve(stand.nobj);
		cmass_harvest_wood_vect.reserve(stand.nobj);
		for(unsigned int i=0;i<stand.nobj;i++) {
			Patch& patch = stand[i];
			cmass_harvest_vect.push_back(check_harvest_cmass(patch));
			cmass_harvest_wood_vect.push_back(check_harvest_cmass(patch, true));
		}
		// Make sure sort_vect is filled
		sort_cmass_harvest();
	}

	int get_index(int rank) {
		return sort_vect[rank].second;
	}
	double get_cmass_harvest(int rank) {
		return cmass_harvest_vect[sort_vect[rank].second];
	}
	double get_cmass_harvest_wood(int rank) {
		return cmass_harvest_wood_vect[sort_vect[rank].second];
	}
};

/// Class for ranking individuals according by diameter etc.
class Rank_individuals {

	Patch& patch;
	std::vector<double> diam_vect;
	std::vector<std::pair<double,int> > sort_vect;

public:

	void sort_data(std::vector<double> v) {

		sort_vect.clear();
		sort_vect.reserve(patch.vegetation.nobj);
		for(unsigned int i=0;i<patch.vegetation.nobj;i++) {
			sort_vect.push_back(std::make_pair(v[i],i));
		}
		std::sort(sort_vect.begin(), sort_vect.end()); 
	}

	void sort_diameter() {
		sort_data(diam_vect);
	}

	Rank_individuals(Patch& patchX) : patch(patchX) {

		diam_vect.reserve(patch.vegetation.nobj);
		for(unsigned int i=0;i<patch.vegetation.nobj;i++) {
			Individual& indiv = patch.vegetation[i];
			diam_vect.push_back(indiv.diam);
		}
		// Make sure sort_vect is filled
		sort_diameter();
	}
	int get_index(int rank) {
		return sort_vect[rank].second;
	}
	double get_diam(int rank) {
		return diam_vect[sort_vect[rank].second];
	}
};

/// Struct for copies of carbon and nitrogen of an individual and associated litter and fluxes resulting from harvest
/// This is needed if we want to harvest only part of a stand, as during land cover change.
struct Harvest_CN {

	double cmass_leaf;
	double cmass_root;
	double cmass_sap;
	double cmass_heart;
	double cmass_debt;
	double cmass_ho;
	double cmass_agpool;
	double cmass_stem;
	double cmass_dead_leaf;
	double debt_excess;
	double nmass_leaf;
	double nmass_root;
	double nmass_sap;
	double nmass_heart;
	double nmass_ho;
	double nmass_agpool;
	double nmass_dead_leaf;
	double nstore_longterm;
	double nstore_labile;
	double max_n_storage;

	double cmass_litter_leaf;
	double cmass_litter_root;
	double cmass_litter_sap;
	double cmass_litter_heart;
	double nmass_litter_leaf;
	double nmass_litter_root;
	double nmass_litter_sap;
	double nmass_litter_heart;
	double acflux_harvest;
	double anflux_harvest;
	double cmass_harvested_products_slow; // May contain original slow product pool value before harvest (if copy_dead_C = true in copy_from_indiv())
	double nmass_harvested_products_slow;

	// The following members are partly overlapping with acflux_harvest; not to be included in balance check functions.
	double acflux_harvest_wood;			// Harvested wood including wood fraction oxidised the same year (1-harvest_slow_frac)
	double acflux_harvest_wood_toprod;	// Always zero before harvest
	double acflux_harvest_tolitter;

	Harvest_CN() {

		cmass_leaf = cmass_root = cmass_sap = cmass_heart = cmass_debt = cmass_ho = cmass_agpool = cmass_stem 
			= cmass_dead_leaf = debt_excess = 0.0;
		nmass_leaf = nmass_root = nmass_sap = nmass_heart = nmass_ho = nmass_agpool = nmass_dead_leaf = nstore_longterm 
			= nstore_labile = max_n_storage = 0.0;
		cmass_litter_leaf = cmass_litter_root = cmass_litter_sap = cmass_litter_heart = 0.0;
		nmass_litter_leaf = nmass_litter_root = nmass_litter_sap = nmass_litter_heart = 0.0;
		acflux_harvest = anflux_harvest = 0.0;
		cmass_harvested_products_slow = nmass_harvested_products_slow = 0.0;
		acflux_harvest_wood = acflux_harvest_wood_toprod = acflux_harvest_tolitter = 0.0;
	}

	/// Copies C and N values from individual and patchpft to struct.
	void copy_from_indiv(Individual& indiv, bool copy_grsC = false, bool copy_dead_C = true) {

		Patch& patch = indiv.vegetation.patch;
		Patchpft& ppft = patch.pft[indiv.pft.id];

		if(copy_grsC) {

			if(indiv.cropindiv) {

				cmass_leaf = indiv.cropindiv->grs_cmass_leaf;
				cmass_root = indiv.cropindiv->grs_cmass_root;

				if(indiv.pft.landcover == CROPLAND) {
					cmass_ho = indiv.cropindiv->grs_cmass_ho;
					cmass_agpool = indiv.cropindiv->grs_cmass_agpool;
					cmass_stem = indiv.cropindiv->grs_cmass_stem;
					cmass_dead_leaf = indiv.cropindiv->grs_cmass_dead_leaf;
				}
			}
		}
		else {

			cmass_leaf = indiv.cmass_leaf;
			cmass_root = indiv.cmass_root;
			cmass_sap = indiv.cmass_sap;
			cmass_heart = indiv.cmass_heart;
			cmass_debt = indiv.cmass_debt;

			if(indiv.pft.landcover == CROPLAND) {
				cmass_ho = indiv.cropindiv->cmass_ho;
				cmass_agpool = indiv.cropindiv->cmass_agpool;
				// cmass_stem = indiv.cropindiv->grs_cmass_stem;			// We can't use grs_cmass here !
				// cmass_dead_leaf = indiv.cropindiv->grs_cmass_dead_leaf;
			}
		}


		nmass_leaf = indiv.nmass_leaf;
		nmass_root = indiv.nmass_root;
		nmass_sap = indiv.nmass_sap;
		nmass_heart = indiv.nmass_heart;
		nstore_longterm = indiv.nstore_longterm;
		nstore_labile = indiv.nstore_labile;
		max_n_storage = indiv.max_n_storage;

		if(indiv.pft.landcover == CROPLAND) {
			nmass_ho = indiv.cropindiv->nmass_ho;
			nmass_agpool = indiv.cropindiv->nmass_agpool;
			nmass_dead_leaf = indiv.cropindiv->nmass_dead_leaf;
		}

		if(copy_dead_C) {

			cmass_litter_leaf = ppft.cmass_litter_leaf;
			cmass_litter_root = ppft.cmass_litter_root;
			cmass_litter_sap = ppft.cmass_litter_sap;
			cmass_litter_heart = ppft.cmass_litter_heart;

			nmass_litter_leaf = ppft.nmass_litter_leaf;
			nmass_litter_root = ppft.nmass_litter_root;
			nmass_litter_sap = ppft.nmass_litter_sap;
			nmass_litter_heart = ppft.nmass_litter_heart;

			// acflux_harvest and anflux_harvest only for output
			cmass_harvested_products_slow = ppft.cmass_harvested_products_slow;
			nmass_harvested_products_slow = ppft.nmass_harvested_products_slow;
		}
	}

	/// Copies C and N values from struct to individual and patchpft living and dead C and N pools. Fluxes are added to patch and patchpft variables.
	/** Use only after a call to copy_from_indiv() before harvest function.
	*/
	void copy_to_indiv(Individual& indiv, bool copy_grsC = false, bool lc_change = false) {

		Patch& patch = indiv.vegetation.patch;
		Patchpft& ppft = patch.pft[indiv.pft.id];

		if(copy_grsC) {

			indiv.cropindiv->grs_cmass_leaf = cmass_leaf;
			indiv.cropindiv->grs_cmass_root = cmass_root;

			if(indiv.pft.landcover == CROPLAND) {
				indiv.cropindiv->grs_cmass_ho = cmass_ho;
				indiv.cropindiv->grs_cmass_agpool = cmass_agpool;
				indiv.cropindiv->grs_cmass_dead_leaf = cmass_dead_leaf;
				indiv.cropindiv->grs_cmass_stem = cmass_stem;
			}
		}
		else {

			indiv.cmass_leaf = cmass_leaf;
			indiv.cmass_root = cmass_root;
			indiv.cmass_sap = cmass_sap;
			indiv.cmass_heart = cmass_heart;
			indiv.cmass_debt = cmass_debt;

			if(indiv.pft.landcover == CROPLAND) {
				indiv.cropindiv->cmass_ho = cmass_ho;
				indiv.cropindiv->cmass_agpool = cmass_agpool;
				// indiv.cropindiv->grs_cmass_dead_leaf = cmass_dead_leaf;	// We can't use grs_cmass here !
				// indiv.cropindiv->grs_cmass_stem = cmass_stem;
			}
		}

		indiv.nmass_leaf = nmass_leaf;
		indiv.nmass_root = nmass_root;
		indiv.nmass_sap = nmass_sap;
		indiv.nmass_heart = nmass_heart;
		indiv.nstore_longterm = nstore_longterm;
		indiv.nstore_labile = nstore_labile;
		indiv.max_n_storage = max_n_storage;

		if(indiv.pft.landcover == CROPLAND) {
			indiv.cropindiv->nmass_ho = nmass_ho;
			indiv.cropindiv->nmass_agpool = nmass_agpool;
			indiv.cropindiv->nmass_dead_leaf = nmass_dead_leaf;
		}

		ppft.cmass_litter_leaf = cmass_litter_leaf;
		ppft.cmass_litter_root = cmass_litter_root;
		ppft.cmass_litter_sap = cmass_litter_sap;
		ppft.cmass_litter_heart = cmass_litter_heart;
		ppft.nmass_litter_leaf = nmass_litter_leaf;
		ppft.nmass_litter_root = nmass_litter_root;
		ppft.nmass_litter_sap = nmass_litter_sap;
		ppft.nmass_litter_heart = nmass_litter_heart;

		if(!lc_change) {
			patch.fluxes.report_flux(Fluxes::HARVESTC, acflux_harvest);	// Put into gridcell.acflux_landuse_change instead at lcc
			patch.fluxes.report_flux(Fluxes::HARVESTN, anflux_harvest);	// Put into gridcell.anflux_landuse_change instead at lcc
		}

		ppft.cmass_wood_harv += acflux_harvest_wood;
		ppft.cmass_wood_harv_toprod += acflux_harvest_wood_toprod;
		ppft.cmass_harv_tolitter += acflux_harvest_tolitter;
	
		// indiv.report_flux(Fluxes::NPP, debt_excess);
		// indiv.report_flux(Fluxes::RA, -debt_excess);

		ppft.cmass_harvested_products_slow = cmass_harvested_products_slow;
		ppft.nmass_harvested_products_slow = nmass_harvested_products_slow;
	}
};

#endif // LPJ_GUESS_MANAGEMENT_H