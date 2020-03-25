///////////////////////////////////////////////////////////////////////////////////////
/// \file soilinput.h
/// \brief Implementation of soil input from text files, soil codes or soil physical properties.
///
///  Created on: 24 nov 2014
/// \author : Stefan Olin
///
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SOILINPUT_H
#define SOILINPUT_H

#include "guess.h"
#include <iterator>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

typedef std::pair<double, double> coord;

/// An input module for soil data.
/** This input module gets soil data from text files, soil codes or soil physical properties.
*/
class SoilInput {
public:
	SoilInput() : 
		searchradius_soil(0),
		STEP(0.5) {
		declare_parameter("searchradius_soil", &searchradius_soil, 0, 100,
			"Search radius for soil input.");
		declare_parameter("searchstep_soil", &STEP, 0, 1,
			"Step size in search function of soil input");
	};

	~SoilInput() {};

	/// Initialises the soil object and decides based on the input file the data source type.
	void init(const char* filename, const std::vector<coord>& gridlist=std::vector<coord>());

	/// Get and set the Soiltype-object in the current Gridcell-object.
	void get_soil(double lon, double lat, Gridcell& gridcell);
	
	double STEP;

private:
	void get_soil_mineral(double lon, double lat, Gridcell& gridcell);

	void get_soil_organic(double lon, double lat, Gridcell& gridcell);

	bool soil_code;
	
	double searchradius_soil;

	coord find_closest_point(double searchradius, coord C);

	struct SoilProperties {

		// empirical parameter in percolation equation (k1) (mm/day)
		double b;
		// volumetric water holding capacity at field capacity minus vol water
		// holding capacity at wilting point (Hmax), as fraction of soil layer
		// depth
		double volumetric_whc_field_capacity;
		// thermal diffusivity (mm2/s) at wilting point (0% WHC)
		double thermal_wilting_point;
		// thermal diffusivity (mm2/s) at 15% WHC
		double thermal_15_whc;
		// thermal diffusivity at field capacity (100% WHC)
		// Thermal diffusivities follow van Duin (1963),
		// Jury et al (1991), Fig 5.11.
		double thermal_field_capacity;
		// wilting point as fraction of depth (calculation method described in
		// Prentice et al 1992)
		double wilting_point;
		// saturation capacity following Cosby (1984)
		double saturation_capacity;
		// sand fraction
		double sand;
		// clay fraction
		double clay;
		// Bulk density
		double bulk_density;
		// pH
		double pH;
		// Organic content (fraction)
		double soil_OC;
		// carbon content (kg C/m2)
 		double soilC;
 		// Porosity
        	double porosity;
	};

	void load_mineral_soils(const char* fname, const std::set<coord>& coords);
	void load_lpj_soilcodes(const char* fname, const std::set<coord>& coords);

	SoilProperties get_lpj(coord c);
	SoilProperties get_mineral(coord c);
	SoilProperties get_lpj_organic_soil();

	struct SoilDataMineral {
		double sand;
		double clay;
		double orgc;
		double bulkdensity;
		double pH;
		double soilC;			// Initial soil C content (kgC/m2), from model output or a database
		double CN;				// C:N ratio from database
	};

	std::map<coord, int> lpj_map;
	std::map<coord, SoilDataMineral> mineral_map;
};

const double data[10][12] = {

	//    0  empirical parameter in percolation equation (k1) (mm/day)
	//    1  volumetric water holding capacity at field capacity minus vol water
	//       holding capacity at wilting point (Hmax), as fraction of soil layer
	//       depth
	//    2  thermal diffusivity (mm2/s) at wilting point (0% WHC)
	//    3  thermal diffusivity (mm2/s) at 15% WHC
	//    4  thermal diffusivity at field capacity (100% WHC)
	//       Thermal diffusivities follow van Duin (1963),
	//       Jury et al (1991), Fig 5.11.
	//    5  wilting point as fraction of depth (calculation method described in
	//       Prentice et al 1992)
	//    6  saturation capacity following Cosby (1984)
	//    7  sand fraction
	//    8  clay fraction
	//    9  volumetric fraction of organic material (m3 m-3) (Hillel, 1998)
	//    10 water held below wilting point, important for heat conductance
	//       From the AGRMET Handbook, 2002 - not used here, use column 5 in data instead
	//    11 porosity from AGRMET Handbook, 2002

	//    0      1      2      3      4      5      6       7       8       9       10      11          soilcode
	//  ------------------------------------------------------------------------------------------
	{ 5.0, 0.110,   0.2, 0.800,   0.4,	0.074,	0.395,	0.90,	0.05,    0.01,    0.029,  0.421 },    // 0	Ice (copied from coarse: Stefan)
	{ 5.0, 0.110,   0.2, 0.800,   0.4,	0.074,	0.395,	0.90,	0.05,    0.01,    0.029,  0.421 },    // 1	Coarse
	{ 4.0, 0.150,   0.2, 0.650,   0.4,	0.184,	0.439,	0.35,	0.15,    0.01,    0.119,  0.464 },    // 2	Medium
	{ 3.0, 0.120,   0.2, 0.500,   0.4,	0.274,	0.454,	0.30,	0.45,    0.01,    0.139,  0.468 },    // 3	Fine
	{ 4.5, 0.130,   0.2, 0.725,   0.4,	0.129,	0.417,	0.60,	0.15,    0.01,    0.047,  0.434 },    // 4	Medium-coarse
	{ 4.0, 0.115,   0.2, 0.650,   0.4,	0.174,	0.425,	0.60,	0.30,    0.01,    0.020,  0.406 },    // 5	Fine-coarse
	{ 3.5, 0.135,   0.2, 0.575,   0.4,	0.229,	0.447,	0.20,	0.30,    0.01,    0.103,  0.465 },    // 6	Fine-medium
	{ 4.0, 0.127,   0.2, 0.650,   0.4,	0.177,	0.430,	0.45,	0.25,    0.01,    0.069,  0.404 },    // 7	Fine-medium-coarse
	{ 9.0, 0.300,   0.1, 0.100,   0.1,	0.200,	0.600,	0.28,	0.12,    0.20,    0.066,  0.800 },    // 8	Organic (values not know for wp), sand and clay are from Parton 2010
	{ 0.2, 0.100,   0.2, 0.500,   0.4,	0.100,	0.250,	0.10,	0.80,    0.01,    0.364,  0.482 }     // 9	Vertisols (values not know for wp)
};

/// Set soil physical properties based in soilcode
void soil_parameters(Soiltype& soiltype, int soilcode);


#endif /* SOILINPUT_H */
