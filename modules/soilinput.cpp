///////////////////////////////////////////////////////////////////////////////////////
/// \file soilinput.cpp
/// \brief Implementation of soil input from text files, soil codes or soil physical properties.
///
///  Created on: 24 nov 2014
/// \author : Stefan Olin
///
/// $Date:  $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "soilinput.h"
#include <fstream>
#include <sstream>


namespace {
bool format_input_header(const char* fname) {
	std::ifstream ifs(fname, std::ifstream::in);
	if (!ifs.good()) {
		fail("SoilInput::init: could not open %s for input", fname);
	}

	std::string line;

	// reads the first line which should be a header
	getline(ifs, line);

	std::istringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> header(begin, end);

	size_t ncols = header.size();
	ifs.close();
	// Return TRUE if the soil data file has three columns (lon, lat and soilcode) and FALSE otherwise (such as when we specify sand, silt, clay, pH etc.) 
	return ncols == 3;
}

}
coord SoilInput::find_closest_point(double searchradius, coord C) {

	// First try the exact coordinate
	bool found = soil_code ? lpj_map.count(C) > 0 :
							 mineral_map.count(C) > 0;
	
	if (found) {
		return C;
	}

	double lon = C.first;
	double lat = C.second;

	if (searchradius == 0) {
		// Don't try to search
		fail("Coordinate at %f, %f could not be found in the soil map, and no search radius specified.\n",lon,lat);
		return C;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(lon * 2) / 2 + 0.25;
	double center_lat = floor(lat * 2) / 2 + 0.25;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy
	// sorting.
	using std::pair;
	using std::make_pair;
	std::vector<pair<double, coord> > search_points;

	
	const double EPS = 1e-15;

	for (double y = center_lon - searchradius; y <= center_lon + searchradius + EPS; y += STEP) {
		for (double x = center_lat - searchradius; x <= center_lat + searchradius + EPS; x += STEP) {
			double xdist = x - lat;
			double ydist = y - lon;
			double dist = sqrt(xdist * xdist + ydist * ydist);

			if (dist <= searchradius + EPS) {
				search_points.push_back(make_pair(dist, make_pair(y, x)));
			}
		}
	}

	// Sort by increasing distance
	std::sort(search_points.begin(), search_points.end());

	// Find closest coordinate which can be found
	for (unsigned int i = 0; i < search_points.size(); i++) {
		coord search_point = search_points[i].second;

		found = soil_code ? lpj_map.count(search_point) > 0 :
							mineral_map.count(search_point) > 0;

		if (found) {
			return search_point;
		}
	}

	fail("Coordinate at %f, %f could not be found in the soil map.\n",lon,lat);
	return C;
}

// Initialises the soil input class and determines the format of the supplied input file.
void SoilInput::init(const char* fname, const std::vector<coord>& gridlist) {

	std::set<coord> coords(gridlist.begin(), gridlist.end());

	soil_code = format_input_header(fname);

	if (soil_code) {
		load_lpj_soilcodes(fname, coords);
	} 
	else {
		load_mineral_soils(fname, coords);
	}
}

// Read in routine for text file with LPJ soil codes.
void SoilInput::load_lpj_soilcodes(const char* fname, const std::set<coord>& coords) {
	std::ifstream ifs(fname, std::ifstream::in);

	std::string line;

	while (getline(ifs, line)) {
		std::istringstream iss(line);

		double lon, lat;
		int classnbr;
		if (iss >> lon >> lat >> classnbr) {
			coord c(lon, lat);
			if (!coords.empty() && coords.count(c) == 0) {
				continue;
			}
			if (classnbr<0 || classnbr>9) {
				fail("SoilInput::init: invalid LPJ soil code (%d) for location"
					" (%g, %g) in file %s", classnbr, lon, lat, fname);
			}
			lpj_map[c] = classnbr;
			if (lpj_map.size() == coords.size()) {
				break;
			}
		}
	}
	ifs.close();
}

// Load mineral soil data from text file with columns describing the data. They need to be named:
// lon, lat, sand, silt, clay, orgc, bulkdensity, ph, soilc, cn.
void SoilInput::load_mineral_soils(const char* fname, const std::set<coord>& coords) {
	std::ifstream ifs(fname, std::ifstream::in);

	std::string line;

	// reads the first line which should be a header
	getline(ifs, line);

	std::transform(line.begin(), line.end(), line.begin(), ::tolower);

	std::istringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> header(begin, end);
	std::vector<std::string>::iterator it = header.begin() + 2;

	// First four variables are intentionally left uninitialised, as those
	// columns are required.
	// int sand_i, clay_i, silt_i, orgc_i, ph_i, bd_i, cn_i, soilc_i = -1;
	int sand_i =-1;
	int clay_i = -1;
	int silt_i = -1; 
	int orgc_i = -1;
	int ph_i = -1;
	int bd_i = -1;
	int cn_i = -1; 
	int soilc_i = -1;

	for (int i = 0; it != header.end(); ++it, ++i) {
		if (*it == "sand") {
			sand_i = i;
		}
		else if (*it == "clay") {
			clay_i = i;
		}	
		else if (*it == "silt") {
			silt_i = i;
		}
		else if (*it == "orgc") {
			orgc_i = i;
		}
		else if (*it == "ph") {
			ph_i = i;
		}
		else if (*it == "bulkdensity") {
			bd_i = i;
		}	
		else if (*it == "cn") {
			cn_i = i;
		}
		else if (*it == "soilc") {
			soilc_i = i;
		}
	}
	if (soilc_i == -1 && iforganicsoilproperties) {
		fail("Error! No Soil C column found in %s. \nTip: do not use iforganicsoilproperties 1 together with a soilmap file without a SoilC column\n", fname);
	}

	// Create a empty vector T with header.size-2 elements
	std::vector<double> T((unsigned int)header.size() - 2);

	while (getline(ifs, line)) {

		std::istringstream iss(line);
		double lon, lat;
		if (iss >> lon >> lat) {	// Get lon and lat from current line.
			coord c(lon, lat);

			if (!coords.empty() && coords.count(c) == 0) {
				continue;
			}

			for (std::vector<double>::iterator it = T.begin(); it != T.end(); ++it) {
				// Add the next values from the line to T
				iss >> *it;
			}

			SoilDataMineral& soildata = mineral_map[c];
			soildata.sand = T[sand_i];
			soildata.clay = T[clay_i];
			soildata.orgc = T[orgc_i];
			soildata.pH = T[ph_i];
			soildata.CN = T[cn_i];
			if (soilc_i != -1) {
				soildata.soilC = T[soilc_i];
			}
			else {
				soildata.soilC = 0.0;
			}
				
			// Not all data sets includes bulk density data, here it is set to a negative number if no column with that name. 
			// TODO, set it to value: 1.6
			if (bd_i<0) {	
				soildata.bulkdensity = (double)bd_i;
			} 
			else {
				soildata.bulkdensity = T[bd_i];
			}

			if (mineral_map.size() == coords.size()) {
				break;
			}
		}
	}

	ifs.close();
}

// Get and set soil properties based on LPJ soil code.
SoilInput::SoilProperties SoilInput::get_lpj(coord c) {

	coord C = find_closest_point(searchradius_soil, c);

	int soilcode = lpj_map[C];
	
	SoilProperties soiltype;
	soiltype.sand = data[soilcode][7];
	soiltype.clay = data[soilcode][8];

	soiltype.b = data[soilcode][0];
	soiltype.volumetric_whc_field_capacity = data[soilcode][1];
	soiltype.thermal_wilting_point = data[soilcode][2];
	soiltype.thermal_15_whc = data[soilcode][3];
	soiltype.thermal_field_capacity = data[soilcode][4];
	soiltype.wilting_point = data[soilcode][5];
	soiltype.saturation_capacity = data[soilcode][6];
	soiltype.pH = 6.5;
	soiltype.soil_OC = 0.05;
	soiltype.soilC = 0.0; // default value (kgC/m2)
	soiltype.porosity = data[soilcode][11];
	return soiltype;
}

// Stores the basic properties of organic soil. 
// Used only if scaling between pure mineral and organic soils is activated.
SoilInput::SoilProperties SoilInput::get_lpj_organic_soil() {

	int soilcode = 8;

	SoilProperties soiltype;
	soiltype.sand = data[soilcode][7];
	soiltype.clay = data[soilcode][8];

	soiltype.b = data[soilcode][0];
	soiltype.volumetric_whc_field_capacity = data[soilcode][1];
	soiltype.thermal_wilting_point = data[soilcode][2];
	soiltype.thermal_15_whc = data[soilcode][3];
	soiltype.thermal_field_capacity = data[soilcode][4];
	soiltype.wilting_point = data[soilcode][5];
	soiltype.saturation_capacity = data[soilcode][6];
	soiltype.pH = 6.5;
	soiltype.soil_OC = 0.05;
	soiltype.soilC = 0.0; // default value (kgC/m2)
	soiltype.porosity = data[soilcode][11];
	return soiltype;
}

// Get and set soil properties based on mineral soil input.
SoilInput::SoilProperties SoilInput::get_mineral(coord c) {
	coord C = find_closest_point(searchradius_soil, c);
	SoilDataMineral& soil = mineral_map[C];
	double silt = 1.0 - soil.sand - soil.clay;

	// Equation 1 from Cosby 1984
	// Psi = Psi_s * (Theta/Theta_s)^b
	// Psi is the pressure head in cm
	// *_s is the values at saturation
	// Theta is the volumetric moisture content in percent
	// Re-arranged to get the Theta
	// Theta = Theta_s * (Psi/Psi_s)^(1/b)

	// from Table 4, Cosby 1984
	double b = 3.10 + 15.7 * soil.clay - 0.3 * soil.sand;

	double logPsi_s = 1.54 - 0.95 * soil.sand + 0.63 * silt;

	// Theta_s in Cosby expressed as %
	double Theta_s = 0.01 * (50.5 - 14.2 * soil.sand - 3.7 * soil.clay);

	double Psi_s = pow(10.0, -logPsi_s);
	double Psi_wilt = pow(10.0, -4.2);
	double Psi_whc = pow(10.0, -2.0);

	double Theta_whc = Theta_s * pow((Psi_whc / Psi_s), 1.0 / b);
	double Theta_wilt = Theta_s * pow((Psi_wilt / Psi_s),1.0 / b);

	// A linear dependence between the percolation coefficient from Haxeltine 1996a
	// and the texture dependent parameter b from Cosby 1984 was established
	// K = 5.87 - 0.29*b

	SoilProperties soiltype;
	soiltype.b = 5.87 - 0.29 * b;

	soiltype.sand = soil.sand;
	soiltype.clay = soil.clay;

	soiltype.volumetric_whc_field_capacity = (Theta_whc - Theta_wilt);
	soiltype.thermal_wilting_point = 0.2;
	soiltype.thermal_15_whc = 0.15 * b + 0.05;
	soiltype.thermal_field_capacity = 0.4;
	soiltype.wilting_point = Theta_wilt;
	soiltype.saturation_capacity = Theta_s;
	soiltype.pH = soil.pH;
	soiltype.soil_OC = soil.orgc;
	soiltype.soilC = soil.soilC;

	if (iforganicsoilproperties) { // Avoid rescaling of porosity twice.
		soiltype.porosity = (1 - soil.orgc) * Theta_s + soil.orgc * organic_porosity; // Eq. 3 Lawrence and Slater, 2008
	}
	else {
		soiltype.porosity = Theta_s;
	}
	
	return soiltype;
}

void SoilInput::get_soil(double lon, double lat, Gridcell& gridcell) {
	
	// Special case if soil properties are presumed to be influenced by the amount of soil carbon. Only valid if not using soil codes.
	if (iforganicsoilproperties) {
		get_soil_organic(lon, lat, gridcell);
	} 
	else {
		get_soil_mineral(lon, lat, gridcell);
	}
}

void SoilInput::get_soil_mineral(double lon, double lat, Gridcell& gridcell) {

	Soiltype& soiltype = gridcell.soiltype;

	coord c(lon, lat);
	SoilProperties soilprop = soil_code ? get_lpj(c) : get_mineral(c);

	soiltype.soilcode = lpj_map[c];
	soiltype.sand_frac = soilprop.sand;
	soiltype.clay_frac = soilprop.clay;
	soiltype.silt_frac = 1.0 - soiltype.sand_frac - soiltype.clay_frac;
	soiltype.perc_base = soilprop.b;
	soiltype.perc_base_evap = soilprop.b;
	soiltype.perc_exp = 2;
	soiltype.gawc[0] = SOILDEPTH_UPPER * soilprop.volumetric_whc_field_capacity;
	soiltype.gawc[1] = SOILDEPTH_LOWER * soilprop.volumetric_whc_field_capacity;
	soiltype.thermdiff_0 = soilprop.thermal_wilting_point;
	soiltype.thermdiff_15 = soilprop.thermal_15_whc;
	soiltype.thermdiff_100 = soilprop.thermal_field_capacity;
	soiltype.gwp[0] = SOILDEPTH_UPPER * soilprop.wilting_point;
	soiltype.gwp[1] = SOILDEPTH_LOWER * soilprop.wilting_point;
	soiltype.gwsats[0] = SOILDEPTH_UPPER * soilprop.saturation_capacity;
	soiltype.gwsats[1] = SOILDEPTH_LOWER *  soilprop.saturation_capacity;
	soiltype.wtot = (soilprop.volumetric_whc_field_capacity + soilprop.wilting_point) * (SOILDEPTH_UPPER + SOILDEPTH_LOWER);

	// Populate the arrays in a more general way
	for (int s = 0; s<NSOILLAYER_UPPER; s++) {
		soiltype.awc[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * soilprop.volumetric_whc_field_capacity;
		soiltype.wp[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * soilprop.wilting_point;
		soiltype.wsats[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * soilprop.saturation_capacity;
	}

	for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {
		soiltype.awc[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * soilprop.volumetric_whc_field_capacity;
		soiltype.wp[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * soilprop.wilting_point;
		soiltype.wsats[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * soilprop.saturation_capacity;
	}

	// Store the soilcode and new soil properties
	soiltype.runon = wetland_runon;
	soiltype.organic_frac = soilprop.soil_OC;
	soiltype.water_below_wp = soilprop.wilting_point;
	soiltype.porosity = soilprop.porosity;
	soiltype.mineral_frac = 1.0 - soiltype.organic_frac - soiltype.porosity;

	if (!ifcentury) {
		// override the default SOM years with 70-80% of the spin-up period
		soiltype.updateSolveSOMvalues(nyear_spinup);
	}
}


void SoilInput::get_soil_organic(double lon, double lat, Gridcell& gridcell) {

	// Use the basic properties of organic soil to update the mineral soil properties, depending on
	// SOC amount and vertical SOC distribution.
	// Scales linearly between pure mineral and organic soils when activated. 

	Soiltype& soiltype = gridcell.soiltype;
	coord c(lon, lat);

	// Get the properties of the mineral soil
	SoilProperties soilpropmineral = soil_code ? get_lpj(c) : get_mineral(c);

	// Determine the soil code, if there is one. If not, set the soilcode to -1.
	int soilcode = soil_code ? lpj_map[c] : -1;

	if (soil_code) { 
		
		// There is no need to update the soil properties below if this already classified as an organic soil type.
		if (soilcode==8) {
			for (int ii = IDX_STD; ii<NLAYERS; ii++) {
				// Save these values for later before updating the thermal properties
				soiltype.org_frac_gridcell[ii - IDX_STD] = 1.0 - organic_porosity;
				soiltype.min_frac_gridcell[ii - IDX_STD] = 0.0;
				soiltype.porosity_gridcell[ii - IDX_STD] = organic_porosity;
			}

			return; 
		}
	}

	// Get the properties of organic soil, i.e. with a soilcode of 8 
	SoilProperties soilproporganic = get_lpj_organic_soil();

	soiltype.sand_frac = soilpropmineral.sand;
	soiltype.clay_frac = soilpropmineral.clay;
	soiltype.silt_frac = 1.0 - soiltype.sand_frac - soiltype.clay_frac;

	//
	// Now update the following properties based on the weighted average of organic properties in the soil.
	//

	// SOC fraction per layer
	double socfrac[NLAYERS];
	// Porosity of each soil layer	
	double porosity[NLAYERS];
	// Organic fraction of the material in each soil layer
	double org_frac[NLAYERS];

	// Total SOC (kgC/m2)
	double soil_tot_soc = 0.0;
	
	// Determine socfrac[NLAYERS] values
	if (!soil_code) {

		soil_tot_soc = soilpropmineral.soilC;

		// Use a beta value to reproduce the global vertical SOC profile shown in Fig 4 of Jobbagy & Jackson (2010)
		const double beta_SOC = 0.976;
		// This value gives SOC fractions as a fraction of the total in the top 100 cm in this Table:
		// Depth range	:	Calculated fraction	:	Jobbagy & Jackson (2000), Fig 4
		// 00 - 20cm	:	0.42				:	0.42
		// 20 - 40cm	:	0.26				:	0.22
		// 40 - 60cm	:	0.16				:	0.16
		// 60 - 80cm	:	0.10				:	0.11
		// 80 - 100cm	:	0.06				:	0.09

		double depth = Dz_soil * CM_PER_MM; // [mm]

		socfrac[IDX_STD] = 1.0 - pow(beta_SOC, depth); // init first layer
		double soc_frac_total = socfrac[IDX_STD];

		for (int ii = IDX_STD+1; ii<NLAYERS; ii++) {
			depth += Dz_soil * CM_PER_MM;
			socfrac[ii] = 1.0 - pow(beta_SOC, depth) - (1.0 - pow(beta_SOC, depth - Dz_soil * CM_PER_MM));
			soc_frac_total += socfrac[ii];
		}

		// Put the residual SOC fraction in lowest soil layer
		socfrac[NLAYERS - 1] += 1.0 - min(soc_frac_total,1.0);
	}

	double carbon_previous_layer = 0.0;
	double soil_c_layer_max = maxSOCdensity * (Dz_soil / MM_PER_M); // kgC m-2 (Dz_soil = 100.0 mm) = 13 kgC m-2

	double mineral_porosity = 0.489 - 0.00126 * 100 * soiltype.sand_frac;	// Eqn 2 (Lawrence & Slater 2008)
	double org_frac_avg = 0.0;
	double org_frac_avg_evap = 0.0;

	double material_min_avg = 0.0; // max = 1.0 - mineral_porosity
	double material_org_avg = 0.0; // max = 0.2

	// Evap layers 
	int num_evaplayers = 2;

	for (int ii = IDX_STD; ii < NLAYERS; ii++) {

		if (!soil_code) {

			double soil_c_layer = soil_tot_soc * socfrac[ii] + carbon_previous_layer; // kgC m-2
			if (soil_c_layer > soil_c_layer_max) {
				carbon_previous_layer = soil_c_layer - soil_c_layer_max;
				soil_c_layer = soil_c_layer_max; // excess carbon saved for the next layer
			}
			else
				carbon_previous_layer = 0.0;

			// Soil C density (kgC/m3) in this soil layer 
			double soil_density_layer = soil_c_layer / (Dz_soil / MM_PER_M); // kgC m-3 (Dz_soil = 100.0 mm)

			// Calculate organic fraction per soil layer using the Lawrence and Slater (2008) approach
			// This is the fraction of the material (i.e. 1 - porosity) that is organic
			org_frac[ii] = soil_density_layer / maxSOCdensity;
			org_frac_avg += org_frac[ii];
			if (ii < IDX_STD + num_evaplayers) {
				org_frac_avg_evap += org_frac[ii];
			}
		}
		else {
			// Use soilpropmineral.soil_OC
			org_frac[ii] = soilpropmineral.soil_OC;
		}

		// Organic fractions > 0 increase the porosity above the mineral value
		porosity[ii] = (1.0 - org_frac[ii]) * mineral_porosity + org_frac[ii] * organic_porosity;

		double material = 1.0 - porosity[ii]; // min value is 0.2

		double material_min = (1.0 - org_frac[ii]) * (1.0 - mineral_porosity); // max = 1.0 - mineral_porosity
		double material_org = org_frac[ii] * (1.0 - organic_porosity); // max = 0.2

		material_min_avg += material_min;
		material_org_avg += material_org;

		// Save these values for later before updating the thermal properties
		soiltype.org_frac_gridcell[ii - IDX_STD] = material_org;
		soiltype.min_frac_gridcell[ii - IDX_STD] = material_min;
		soiltype.porosity_gridcell[ii - IDX_STD] = porosity[ii];
	}

	if (soil_code) {
		org_frac_avg = soilpropmineral.soil_OC;
		org_frac_avg_evap = soilpropmineral.soil_OC;
	}
	else {
		org_frac_avg /= (NLAYERS - IDX_STD);
		org_frac_avg_evap /= (double)num_evaplayers;
	}

	double min_frac_avg = 1.0 - org_frac_avg;
	double min_frac_avg_evap = 1.0 - org_frac_avg_evap;

	material_min_avg /= (NLAYERS-IDX_STD);
	material_org_avg /= (NLAYERS-IDX_STD);

	// First calculate the organic fractions in the top and the bottom layers
	double org_frac_upper = 0.0;
	for (int s = 0; s < NSOILLAYER_UPPER; s++) {
		org_frac_upper += org_frac[IDX_STD + s] / NSOILLAYER_UPPER;
	}
	double min_frac_upper = 1.0 - org_frac_upper;

	double org_frac_lower = 0.0;
	for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {
		org_frac_lower += org_frac[IDX_STD + s] / NSOILLAYER_LOWER;
	}
	double min_frac_lower = 1.0 - org_frac_lower;

	soiltype.perc_base = org_frac_avg * soilproporganic.b + min_frac_avg * soilpropmineral.b;
	soiltype.perc_base_evap = org_frac_avg_evap * soilproporganic.b + min_frac_avg_evap * soilpropmineral.b;
	soiltype.perc_exp = 2;

	soiltype.gawc[0] = SOILDEPTH_UPPER * (min_frac_upper * soilpropmineral.volumetric_whc_field_capacity + 
		org_frac_upper * soilproporganic.volumetric_whc_field_capacity);
	soiltype.gawc[1] = SOILDEPTH_LOWER * (min_frac_lower * soilpropmineral.volumetric_whc_field_capacity + 
		org_frac_lower * soilproporganic.volumetric_whc_field_capacity);
	
	soiltype.thermdiff_0 = org_frac_avg * soilproporganic.thermal_wilting_point + min_frac_avg * soilpropmineral.thermal_wilting_point;
	soiltype.thermdiff_15 = org_frac_avg * soilproporganic.thermal_15_whc + min_frac_avg* soilpropmineral.thermal_15_whc;
	soiltype.thermdiff_100 = org_frac_avg * soilproporganic.thermal_field_capacity + min_frac_avg * soilpropmineral.thermal_field_capacity;

	soiltype.gwp[0] = SOILDEPTH_UPPER * (min_frac_upper * soilpropmineral.wilting_point + 
		org_frac_upper * soilproporganic.wilting_point);
	
	soiltype.gwp[1] = SOILDEPTH_LOWER * (min_frac_lower * soilpropmineral.wilting_point + 
		org_frac_lower * soilproporganic.wilting_point);

	soiltype.gwsats[0] = SOILDEPTH_UPPER * (min_frac_upper * soilpropmineral.saturation_capacity + 
		org_frac_upper * soilproporganic.saturation_capacity);

	soiltype.gwsats[1] = SOILDEPTH_LOWER * (min_frac_lower * soilpropmineral.saturation_capacity + 
		org_frac_lower * soilproporganic.saturation_capacity);

	soiltype.wtot = soiltype.gawc[0] + soiltype.gawc[1] + soiltype.gwp[0] + soiltype.gwp[1];

	// Populate the arrays in a more general way
	for (int s = 0; s < NSOILLAYER_UPPER; s++) {
		soiltype.awc[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * ((1.0 - org_frac[IDX_STD + s]) * soilpropmineral.volumetric_whc_field_capacity + 
			org_frac[IDX_STD + s] * soilproporganic.volumetric_whc_field_capacity);

		soiltype.wp[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * ((1.0 - org_frac[IDX_STD + s]) * soilpropmineral.wilting_point + 
			org_frac[IDX_STD + s] * soilproporganic.wilting_point);

		soiltype.wsats[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * ((1.0 - org_frac[IDX_STD + s]) * soilpropmineral.saturation_capacity + 
			org_frac[IDX_STD + s] * soilproporganic.saturation_capacity);
	}

	for (int s = NSOILLAYER_UPPER; s < NSOILLAYER; s++) {
		soiltype.awc[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * ((1.0 - org_frac[IDX_STD + s]) * soilpropmineral.volumetric_whc_field_capacity + 
			org_frac[IDX_STD + s] * soilproporganic.volumetric_whc_field_capacity);

		soiltype.wp[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * ((1.0 - org_frac[IDX_STD + s]) * soilpropmineral.wilting_point + 
			org_frac[IDX_STD + s] * soilproporganic.wilting_point);

		soiltype.wsats[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * ((1.0 - org_frac[IDX_STD + s]) * soilpropmineral.saturation_capacity + 
			org_frac[IDX_STD + s] * soilproporganic.saturation_capacity);
	}

	// Store the soilcode and new soil properties
	soiltype.runon = wetland_runon;
	soiltype.soilcode = soilcode -1;								// -1 if soil_code is false
	soiltype.water_below_wp = soilpropmineral.wilting_point;
	// These values are not used when we use the organic fraction to determine soil properties
	// However, here we update them with the average values
	soiltype.organic_frac = material_org_avg;
	soiltype.mineral_frac = material_min_avg;
	soiltype.porosity = 1.0 - soiltype.organic_frac - soiltype.mineral_frac;

	if (!ifcentury) {
		// override the default SOM years with 70-80% of the spin-up period
		soiltype.updateSolveSOMvalues(nyear_spinup);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// SOILPARAMETERS
// May be called from input/output module to initialise stand Soiltype objects when
// soil data supplied as LPJ soil code rather than soil physical parameter values
void soil_parameters(Soiltype& soiltype, int soilcode) {

	// DESCRIPTION
	// Derivation of soil physical parameters given LPJ soil code

	// INPUT AND OUTPUT PARAMETER
	// soil = patch soil

	// exponent in percolation equation [k2; LPJF]
	// (Eqn 31, Haxeltine & Prentice 1996)
	// Changed from 4 to 2 (Sitch, Thonicke, pers comm 26/11/01)
	const double PERC_EXP = 2.0;

	if (soilcode<0 || soilcode>9) {
		fail("soil_parameters: invalid LPJ soil code (%d)", soilcode);
	}

	if (textured_soil) {
		soiltype.sand_frac = data[soilcode][7];
		soiltype.clay_frac = data[soilcode][8];
	}
	else {
		// Using fixed values from Parton et al. (2010)
		soiltype.sand_frac = 0.28;
		soiltype.clay_frac = 0.12;
	}

	soiltype.silt_frac = 1.0 - soiltype.sand_frac - soiltype.clay_frac;
	soiltype.perc_base = data[soilcode][0];
	soiltype.perc_base_evap = data[soilcode][0];
	soiltype.perc_exp = PERC_EXP;
	soiltype.gawc[0] = SOILDEPTH_UPPER * data[soilcode][1];
	soiltype.gawc[1] = SOILDEPTH_LOWER * data[soilcode][1];
	soiltype.thermdiff_0 = data[soilcode][2];
	soiltype.thermdiff_15 = data[soilcode][3];
	soiltype.thermdiff_100 = data[soilcode][4];
	soiltype.gwp[0] = SOILDEPTH_UPPER * data[soilcode][5];
	soiltype.gwp[1] = SOILDEPTH_LOWER * data[soilcode][5];
	soiltype.gwsats[0] = SOILDEPTH_UPPER * data[soilcode][6];
	soiltype.gwsats[1] = SOILDEPTH_LOWER * data[soilcode][6];
	soiltype.wtot = (data[soilcode - 1][1] + data[soilcode][5]) * (SOILDEPTH_UPPER + SOILDEPTH_LOWER);

	// Populate the arrays in a more general way
	for (int s = 0; s < NSOILLAYER_UPPER; s++) {
		soiltype.awc[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * data[soilcode][1];
		soiltype.wp[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * data[soilcode][5];
		soiltype.wsats[s] = SOILDEPTH_UPPER / NSOILLAYER_UPPER * data[soilcode][6];
	}

	for (int s = NSOILLAYER_UPPER; s < NSOILLAYER; s++) {
		soiltype.awc[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * data[soilcode][1];
		soiltype.wp[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * data[soilcode][5];
		soiltype.wsats[s] = SOILDEPTH_LOWER / NSOILLAYER_LOWER * data[soilcode][6];
	}

	// Store the soilcode and new soil properties
	soiltype.runon = wetland_runon;
	soiltype.soilcode = soilcode;
	soiltype.organic_frac = data[soilcode][9];
	soiltype.water_below_wp = data[soilcode][5];
	soiltype.porosity = data[soilcode][11];
	soiltype.mineral_frac = 1.0 - soiltype.organic_frac - soiltype.porosity;

	if (!ifcentury) {
		// override the default SOM years with 70-80% of the spin-up period
		soiltype.updateSolveSOMvalues(nyear_spinup);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Jobbagy, E. G. & Jackson, R. B. (2000) The vertical distribution of soil carbon
//  and its relation to climate and vegetation. Ecological Applications 10(2): 423-436
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical exploration
//   of the relationships of soil moisture characteristic to the physical properties of soil.
//   Water Resources Research, 20: 682-690.
