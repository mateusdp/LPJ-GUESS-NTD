///////////////////////////////////////////////////////////////////////////////////////
/// \file soil.h
/// \brief Constants and parameters used in Arctic and wetland code, with references.
/// NB: The class Soil and its member functions and variables are declared in guess.h, 
/// while its member functions are implemented in soil.cpp and in soilmethane.cpp.
/// 
/// \author Paul Miller
/// $Date: 2019-03-10 16:34:14 +0100 (Sun, 10 Mar 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_SOIL_H
#define LPJ_GUESS_SOIL_H

// DEBUGGING BOOLEANS

const bool DEBUG_SOIL_WATER = false;
const bool DEBUG_SOIL_TEMPERATURE = false;
const bool DEBUG_METHANE = false;


// CONSTANTS

/// number of total soil layers. Must be at least NSOILLAYER + NLAYERS_SNOW + 1
const int NLAYERS = 22;

/// Number of soil layers for soil temperature/water calculations. Typically 15 10mm layers, making up the 150cm-deep soil column
const int NSOILLAYER = 15;	// rootdist in .ins file must have NSOILLAYER components

/// Number of soil layers used in LPJ-GUESS v4.0 for soil water calculations. Typically 2 layers, 500mm + 10000mm, making up the 150cm-deep soil column
const int NSOILLAYER_SIMPLE = 2;

/// number of depths at which we want the soil T output in outannual
const int SOILTEMPOUT = NSOILLAYER;

/// index of the first soil layer
const int IDX_STD = NLAYERS-NSOILLAYER;

/// number of padding layers in the soil
const int PAD_LAYERS = 5;

/// number of total soil layers in the acrotelm
const int NACROTELM = 3;

/// number of total soil layers in the catotelm
const int NCATOTELM = 12;

/// Number of sublayers in the acrotelm, i.e. 1cm layers with this = 30
const int NSUBLAYERS_ACRO = 30;

/// Maximum number of layers - for soil temperature calculations
const int active_layersmax = NLAYERS + PAD_LAYERS;

/// Total depth [mm] of padding layers - set to 8000 when analyticalSolutionTest = true
const double PAD_DEPTH = 48000;

/// Depth of each mineral soil layer [mm]
const double Dz_soil = 100.0;

/// Depth of each acrotelm soil layer [mm]
const double Dz_acro = 100.0;

/// Depth of each catotelm soil layer [mm]
const double Dz_cato = 100.0;

/// Max height of standing water [mm]
const double maxh = 0.0;

/// Slope of soil profile
const double soil_slope = -.37;

/// Maximum density of soil organic carbon [KgC/m3]
const double maxSOCdensity = 130.0; // From Lawrence and Slater, 2008

// SNOW PARAMETERS

/// snow density at start of the snow season [kg m-3] (Wania et al. (2009a) have 150, Best et al. (2012) have 50)
const double snowdens_start = 275.0;

/// snow density at end of the snow season [kg m-3] (Wania et al. 2009a)
const double snowdens_end = 500.0;

/// maximum number of snow layers allowed (<= 5)
const int NLAYERS_SNOW = 5;

/// ice density [kg m-3] - CLM value
const double ice_density = 917.0;

/// water density [kg m-3]
const double water_density = 1000.0;

// POROSITIES

/// Porosity of organic material
const double organic_porosity = 0.8; // From Lawrence and Slater (2008) have 0.9, but 0.8 is consistent with organic soil code 8

/// catotelm porosity
const double catotelm_por = 0.92;

/// acrotelm porosity
const double acrotelm_por = 0.98;

/// Gas fraction in peat
const double Fgas = 0.00; // Possible to reintroduce - was 0.08 in Wania et al (2010)

/// Wilting point in peat
const double peat_wp = 0.066;

/// First year when phase change is allowed
const int FIRST_FREEZE_YEAR = 90;

/// time step [day]
const double Dt = 1;

// HEAT CAPACITIES

/// heat capacity of air [J m-3 K-1] - Bonan (2002)
const double Cp_air = 1200;

/// heat capacity of water [J m-3 K-1] - Bonan (2002)
const double Cp_water = 4180000;
 
/// heat capacity of ice [J m-3 K-1] = 2117.27 [J kg-1 K-1] * 917 [kg m-3 (ice density)] - CLM
const double Cp_ice = 1941537;		 

/// heat capacity of organic matter [J m-3 K-1]
const double Cp_org = 2500000;

/// heat capacity of dry peat (J m-3 K-1) - Bonan (2002), 0% water
const double Cp_peat = 580000;

/// heat capacity of mineral soil [J m-3 K-1] - Bonan (2002)
const double Cp_min = 2380000;

/// heat capacity of moss [J/ m-3 K-1] Ekici et al. (2015)
const double Cp_moss = 2500000;

// NOTE: using the Cp_org values for Cp_peat and Korg for Kpeat does not
// seem to influence the upper soil layer Ts, but increases the range
// in the lower layers (2m)

// THERMAL CONDUCTIVITIES
// Values are from Hillel (1982) unless otherwise stated

/// thermal conductivity of air [W m-1 K-1]
const double Kair = 0.025;

/// thermal conductivity of water [W m-1 K-1]
const double Kwater = 0.57;

/// thermal conductivity of ice [W m-1 K-1]
const double Kice = 2.2;

/// thermal conductivity of organic matter [W m-1 K-1]
const double Korg = 0.25;

/// thermal conductivity of dry peat [W m-1 K-1] - Bonan (2002)
const double Kpeat = 0.06;

/// thermal conductivity of mineral soil [W m-1 K-1] - Wania et al. (2009a)	 
const double Kmin = 2.0;

/// thermal conductivity of moss [W m-1 K-1]
const double Kmoss = 0.25;

/// latent heat of fusion (Granberg et al. 1999) [J m-3]
const double Lheat = 3.34E8;	

/// A FILL-IN for layers above layer0.
const double MISSING_VALUE = -9999.0;      

/// The number of subdaily timestep loops to perform. 
// Numerical instabilities can arise when there are thinner layers of snow and litter, so this should be > 1. 
const int TIMESTEPS = 2;

/// Thickness of the topmost air layer [m]
// Note - keep air thickness low when using cnstep_full
const double AIR_THICKNESS = 100.0;


////////////////////////////////////////////////////////////////////////

// SOIL CONSTANTS 

/// reduce infiltration and percolation rates sharply when there is ice present in the soil
// See CLM4.5 - Swenson et al. (2012)
const bool ICE_IMPEDANCE = false;
	
/// min temp [deg C] for heterotrophic decomposition.
// Clein & Schimel (1995) use -4 degC
const double MIN_DECOMP_TEMP = -8.0;		
	
/// max CO2 [mimol/L] available to mosses in the acrotelm - see Wania et al (2009b)
// value taken from Smolders et al. 2001 which give an average CO2 conc. of 70 sites as 934 mimol L-1
const double PORE_WATER_CO2 = 934.0; 

// METHANE CONSTANTS 

/// optimised moisture response under inundation (Wania et al. 2010, Table 5)
const double RMOIST = 0.4;

/// Frolking et al (2001, 2010), Ise et al. (2008)
const double RMOIST_ANAEROBIC=0.025; 
	
/// time step for gas diffusion calculations [day]
const double Dt_gas = 0.01;
	
/// CH4:CO2 ratio for peatland soils (> PEATLAND_WETLAND_LATITUDE_LIMIT N) - See Wania et al (2010)
// Wania et al. (2010) optimal value: 0.1 (see Table 4). 
// McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.25, after optimisation
const double CH4toCO2_peat = 0.085; 
	
/// CH4:CO2 ratio for inundated soils (< PEATLAND_WETLAND_LATITUDE_LIMIT N) - See Spahni et al. (2011)
//const double CH4toCO2_inundated = 0.024; // SC1 value in Spahni et al. SC2 is 0.0415
const double CH4toCO2_inundated = 0.027; // Updated from Spahni et al. (2011) to match global emissions

/// density of water [kg m-3]
const double rho_H2O = 1000.0;              
	
/// acceleration due to gravity [m s-2]
const double gravity = 9.81;

/// Molecular mass of CH4 [g mol-1]
const double mr_CH4 = 16.0;

/// Molecular mass of carbon [g mol-1]
const double mr_C = 12.0;

/// molecular weight of water
const double mr_h2o = 18.0;

/// universal gas constant [J mol-1 K-1]
const double R_gas = 8.314472;

/// standard atmospheric pressure [Pa]
const double atm_press = 101325.0;

/// coefficient for the calculation of the gas transport velocity, given in Riera et al. 1999
const double n_coeff = -0.5;

/// wind speed at 10m height [m s-1]
const double U10 = 0.0;
      
///  Henry's Law constants [L atm mol-1] at 298.15K. Wania et al. (2010), Table 2
const double henry_k_CO2 = 29.41;
const double henry_k_CH4 = 714.29;
const double henry_k_O2 = 769.23;
     
/// Constants [K] for CO2, CH4 and O2 for calculation of Henry's coefficient cited by Sander (1999). Wania et al. (2010), Table 2
const double henry_C_CO2 = 2400.0;
const double henry_C_CH4 = 1600.0;
const double henry_C_O2 = 1500.0;

/// partial pressure of CH4 above water
const double pp_CH4 = 1.7; // micro atm

/// partial pressure of O2 above water (value consistent with PO2 in canexch.h)
const double pp_O2 = 209000; // micro atm

/// when ebullition occurs, the volumetric gas content (VGC) will drop to this level [unitless]
const double vgc_low = 0.145;

///	ebullition occurs, when the volumetric gas content (VGC) exceeds this level [unitless, m3/m3]
const double vgc_high = 0.15;          

/// CH4 fraction of gas bubbles [unitless]
const double bubble_CH4_frac = 0.57;

/// Fraction of oxygen used to oxidise CH4
// Wania et al. (2010) optimal value: 0.5 (see Table 5). 
// McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.9, after optimisation
const double oxid_frac = 0.5;

/// Fraction of ANPP used to calculate number of tillers
const double ag_frac = 0.4;

/// Radius of an average tiller [m]
// (tiller_radius = 0.004)  ! Schimel (1995) - Average over E. angustifolium (diam=7.9mm)
// and C. aquatilis (diam=3.8mm)
// Wania et al. (2010) optimal value: 0.003mm (see Table 5). 
// McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.0035, after optimisation
const double tiller_radius = 0.0035;       


/// Tiller porosity
// (tiller_por = 0.6) ! Wetland plants book, eds. Cronk and Fennessy, p.90, values for 2 Erioph. spp.
const double tiller_por = 0.7; // Wania et al. (2010) optimal value: 0.7 (see Table 5)
   
/// C content of biomass
const double c_content = 0.45;

/// atomic mass of carbon [g/mol]
const double atomiccmass = 12.0;

/// Individual tiller weight [g C]
const double tiller_weight = 0.22;
     
/// a threshold factor for a minimum water content in the layer [unitless]
const double water_min = 0.1;

/// Latitude (N). North of this and PEATLAND stands are treated as peatland as in Wania et al. (2009a, 2009b, 2010)
// But south of this, then the PEATLAND stands are irrigated to avoid water stress, and can be a source of methane
const double PEATLAND_WETLAND_LATITUDE_LIMIT = 40.0;

// INLINE FUNCTIONS

inline void tridiag(int n, long double a[], long double b[], long double c[], long double r[], long double u[]) {

	// Tridiagonal system solver from Numerical Recipes.

	double gam[active_layersmax];
	double bet;

	bet = b[0];

	u[0] = r[0] / bet;

	for (int j = 1; j<n; j++) {
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];

		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}

	for (int j = (n - 2); j >= 0; j--) {
		u[j] -= gam[j + 1] * u[j + 1];
	}
}


///////////////////////////////////////////////////////////////////////////////

// References

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

#endif //LPJ_GUESS_SOIL_H
