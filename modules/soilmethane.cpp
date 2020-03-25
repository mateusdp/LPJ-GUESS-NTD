///////////////////////////////////////////////////////////////////////////////////////
/// \file methane.cpp
/// \brief Methane calculations and transport.
/// Implementation of member functions of class Soil.
/// The class Soil and its member functions and variables are declared in guess.h
///
/// \author Paul Miller
/// $Date: 2017-05-203 13:32:28 +0200 (Wed, 03 May 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include "soilmethane.h"
#include "guess.h"
#include "soil.h"


//////////////////////////////////////////////////////////////////////////////////////////////////
// METHANE
//////////////////////////////////////////////////////////////////////////////////////////////////


/// Calculate methane emissions and update methane concentrations in soil layers. 
/**
	Emissions if and only if ifmethane true (in global.ins) and this is a PEATLAND stand
	and we have passed the initial years of the spinup.

	Called daily from framework.

*/
void methane_dynamics(Patch& patch) {

	if (ifmethane && patch.stand.landcover == PEATLAND && date.year > nyear_spinup-100) {
		patch.soil.methane(true);
	} 
	else {
		patch.soil.methane(false);
	}
}

/// Crank-Nicholson timestepper algorithm for gas diffusion equation.
void cnstepgas(int layer0, double Di[NLAYERS], double dz[NLAYERS], double surf_conc, double dt, 
			double conc[NLAYERS]) {

	// UNITS:
	// Di			m2 d-1
	// dz			m
	// surf_conc	mmol m-3
	// dt			d
	// conc			mmol m-3

	// A fill-in for layers above layer0.
	double MISSING_VALUE;

	// Layer counters: note that there are two different layer counting
	// schemes used, one for the input and output parameters (vectors of
	// length NLAYERS) and one for the values used in the Crank-Nicholson
	// solver (vectors of length active_layers)

	int layer, lidx;
	int active_layers;
	double dplus;
	double dminus;
	double dz_factor;
	double Cplus; 
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double cohere, cominus, coplus;

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
	dplus=0.0;
	dminus=0.0;
	dz_factor=0.0;
	Cplus=0.0;
	Cp_minus=0.0;
	dzhere=0.0;
	dzminus=0.0;
	dzplus=0.0;
	cohere=0.0;
	cominus=0.0;
	coplus=0.0;

	for (int l = 0; l < active_layersmax; l++) {
		diag[l] = 0.0;
		left[l] = 0.0;
		right[l] = 0.0;
		rhs[l] = 0.0;
		solution[l] = 0.0;
	}

    // --- CODE STARTS HERE ---

	MISSING_VALUE = -9999.0;

	active_layers = NLAYERS - layer0; //layer0 = IDX

	//   BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// End members for off-diagonal elements.
	left[0] = 0.0;
	right[active_layers-1] = 0.0;

	// Process the active layers. 
	// The first time (lidx and layer = layer0) corresponds to the surface layer

	for (lidx = 1; lidx<=active_layers; lidx++) {

		// Deal with different layer counting schemes.
		layer = lidx + layer0 - 1; 
		// Minimum is layer0
		// Maximum is NLAYERS - 1

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.  

		// D+

		if (layer==NLAYERS -1) 
			dplus = 0.0; // BC2 - Bottom layer diffusion clamped to 0
		else 
			dplus = 0.5 * (Di[layer] + Di[layer + 1]);

		// D-

		if (layer==layer0) 
			dminus = Di[layer]; // top layer
		else
			dminus = 0.5 * (Di[layer] + Di[layer - 1]); // soil layers	
 
		// --- HERE ---
		if (layer < NLAYERS) {
			// all soil layers
			dzhere = dz[layer];
			cohere = conc[layer];
		}

		// --- PLUS ---
		if (layer < NLAYERS - 1) {
			// all soil layers apart from the bottom soil layer
			dzplus = dz[layer + 1];
			coplus = conc[layer + 1];
		}

		// --- MINUS ---
		if (layer==layer0) {
			// top layer
			dzminus = dz[layer];
			cominus = conc[layer];
		}
		else if (layer < NLAYERS) {
			// all soil layers and the first padding layer
			dzminus = dz[layer - 1];
			cominus = conc[layer - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.
 
		// DIAG
		if (lidx==1) 
			diag[0] = 1.0; // BC1 - top layer should be (1,0,...,0)
		else
			diag[lidx-1] = 1.0 + Cplus + Cp_minus;
 
		// LEFT & RIGHT
		if (lidx < active_layers) {

			if (lidx > 1)
				right[lidx-1] = -Cplus;
			else
				right[lidx-1] = 0.0; // i.e. BC1, where top layer == (1,0,..,0)

			// left[0] is set above.
			if (lidx > 1)
				left[lidx-1] = -Cp_minus;			
		}

		if (lidx == active_layers) 
			left[lidx-1] = -Cp_minus;
		 
		// RHS
		// Calculate right hand side vector values.
		if (lidx==1) 
			rhs[0] = surf_conc;
		else if (lidx==active_layers) // Cplus == 0 here anyway
			rhs[lidx-1] = (1.0 - Cp_minus) * cohere + Cp_minus * cominus;
		else
			rhs[lidx-1] = (1.0 - Cplus - Cp_minus) * cohere +
			Cplus * coplus + Cp_minus * cominus;
	} // end for

	//   SOLVE TRIDIAGONAL SYSTEM
	tridiag(active_layers, left, diag, right, rhs, solution);

	// Numerical tests

	if (DEBUG_METHANE) {

		// Test 1:
		int testrow = 10; 
		double checksum = left[testrow]*solution[testrow-1] + diag[testrow]*solution[testrow] +
		right[testrow]*solution[testrow+1] - rhs[testrow]; 

		const double MAX_ERR = 0.000001;
	
		if (fabs(checksum) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s\n","Bad checksum after cnstepgas - tridiag - test1");

		// Test 2:
		// Every entry in rowsum should == 1.0 
		long double rowsum[active_layersmax];
		checksum = 0.0;	  
		for (lidx = 0; lidx<active_layers; lidx++) {
			rowsum[lidx] = left[lidx] + diag[lidx] + right[lidx];
			checksum += rowsum[lidx] / double(active_layers);
		}

		if (fabs(checksum-1.0) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s\n","Bad checksum after cnstepgas - tridiag - test2");
	}

	//   FORMAT OUTPUT PARAMETERS

	// Transfer the solution to the concentration array.
	for (lidx = 0; lidx<NLAYERS; lidx++) {

		if (lidx < layer0) 
			conc[lidx] = MISSING_VALUE;
		else
			conc[lidx] = solution[lidx - layer0];
	}
}


bool Soil::calculate_gas_diffusivities(double dCH4[NLAYERS], double dCO2[NLAYERS], double dO2[NLAYERS]) {

	// Calculate and return diffusivities, in units of m2 d-1
	// Called each day
	// See Wania et al. (2010) - Sec 2.5

	double D_CH4_water[NLAYERS];
	double D_CO2_water[NLAYERS];
	double D_O2_water[NLAYERS];

	double D_CH4_air[NLAYERS];
	double D_CO2_air[NLAYERS];
	double D_O2_air[NLAYERS];

	for (int ii=0; ii<NLAYERS; ii++) {
		D_CH4_water[ii] = 0.0;
		D_CH4_air[ii] = 0.0;
		dCH4[ii] = 0.0;
		D_CO2_water[ii] = 0.0;
		D_CO2_air[ii] = 0.0;
		dCO2[ii] = 0.0;
		D_O2_water[ii] = 0.0;
		D_O2_air[ii] = 0.0;
		dO2[ii] = 0.0;
	}

	const double scale = 0.00001;

	for (int ii=IDX; ii<NLAYERS; ii++) {

		// AIR - UNITS: 10-4 m2 s-1 == cm2 s-1
		// Lerman (1979)

		double layerT = T_soil[ii];

		// Wania et al. (2010), Eqn. 10
		D_CH4_air[ii] = 0.1875 + 0.00130 * layerT; 
		D_CO2_air[ii] = 0.1325 + 0.00090 * layerT;
		D_O2_air[ii]  = 0.1759 + 0.00117 * layerT;

		// WATER - UNITS: 10-4 m2 s-1 == cm2 s-1
		// Wania et al. (2010), Eqn. 9 - where the diffusivities are in units of 10**-9 m2 s-1. 
		// * by scale to get 10-4 m2 s-1
		D_CH4_water[ii] = (0.9798 + 0.02986 * layerT + 0.0004381 * layerT * layerT) * scale;
		D_CO2_water[ii] = (0.9390 + 0.02671 * layerT + 0.0004095 * layerT * layerT) * scale;
		D_O2_water[ii]  = (1.1720 + 0.03443 * layerT + 0.0005048 * layerT * layerT) * scale;
	}

	// ACROTELM diffusivities 
	for (int ii = IDX; ii < IDX + NACROTELM; ii++) {

		if (Frac_air[ii] > 0.05) {
			// Wania et al. (2010), Eqn. 11 
			double airpow = pow(Frac_air[ii],10/3) / pow(acrotelm_por,2.0);
			dCH4[ii] = airpow * D_CH4_air[ii];
			dCO2[ii] = airpow * D_CO2_air[ii];
			dO2[ii]  = airpow * D_O2_air[ii];				
		} 
		else {
			// Water diffusivities dominate
			dCH4[ii] = D_CH4_water[ii];
			dCO2[ii] = D_CO2_water[ii];
			dO2[ii] = D_O2_water[ii];			
		}
	}

	// CATOTELM - assumed to be always saturated, so use water diffusivities
	for (int ii = IDX + NACROTELM; ii < NLAYERS; ii++) {
		dCH4[ii] = D_CH4_water[ii];
		dCO2[ii] = D_CO2_water[ii];
		dO2[ii] = D_O2_water[ii];			
	}

	for (int ii=IDX; ii<NLAYERS; ii++) {

		// Convert from cm2 s-1 to m2 d-1
		dCH4[ii] *= SECS_PER_DAY / CM2_PER_M2;
		dCO2[ii] *= SECS_PER_DAY / CM2_PER_M2;
		dO2[ii] *= SECS_PER_DAY / CM2_PER_M2;
	}

	return true;
}


bool Soil::update_daily_gas_parameters() {

	// Update gas transport velocities [m d-1] and equilibrium gas concentrations [mmol m-3]
	// for CH4, CO2 and O2
	// Called daily by Soil::methane()
	// See Wania et al. (2010) - Sec 2.5

	double surfT = T_soil[IDX];

	/* GAS TRANSPORT VELOCITIES */

	// gas transport velocity of SF6.
	// Wania et al. (2010) - Eqn. 6
	double k_600 = 2.07 + 0.215 * pow(U10,1.7); // cm h-1

	// Schmidt number of O2
	// Wania et al. (2010) - Eqn. 7
	double ScO2 = 1800.6 - 120.1 * surfT + 3.7818 * pow(surfT,2) - 0.047608 * pow(surfT,3);
	
	// gas transport velocity of O2 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	k_O2 = k_600 * pow(ScO2/600.0, n_coeff);
	k_O2 *= 24/CM_PER_M; // cm h-1 to m d-1
	
	// Schmidt number of CH4
	// Wania et al. (2010) - Eqn. 7
	double ScCH4 = 1898.0 - 110.1 * surfT + 2.834 * pow(surfT,2) - 0.02791 * pow(surfT,3);
	
	// gas transport velocity of CH4 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	k_CH4 = k_600 * pow(ScCH4/600.0, n_coeff);
	k_CH4 *= 24/ CM_PER_M; // cm h-1 to m d-1

	// Schmidt number of CH4
	// Wania et al. (2010) - Eqn. 7
	double ScCO2 = 1911.0 - 113.7 * surfT + 2.967 * pow(surfT,2) - 0.02943 * pow(surfT,3);
	
	// gas transport velocity of CO2 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	k_CO2 = k_600 * pow(ScCO2/600.0, n_coeff);
	k_CO2 *= 24/CM_PER_M; // cm h-1 to m d-1

	/* EQUILIBRIUM GAS CONCENTRATIONS */
	// See Wania et al. (2010) - Eqn. 8

	double deg25 = K2degC + 25.0; // 25 degrees C

	double henry_coeff_O2 = henry_k_O2 * exp(-1.0 * henry_C_O2 * (1.0 / (surfT + K2degC) - 1.0 / deg25));
	// [L atm mol-1]

	// pp_gas/MM2_PER_M2 converts to atm units

	Ceq_O2 = pp_O2 / MM2_PER_M2 / henry_coeff_O2; // mol L-1
	Ceq_O2 *= MM2_PER_M2; // mol L-1 to mmol m-3

	double henry_coeff_CO2 = henry_k_CO2 * exp(-1.0 * henry_C_CO2 * (1.0 / (surfT + K2degC) - 1.0 / deg25));
	// [L atm mol-1]

	// Use this gridcell's CO2 concentration, not a fixed value as in Wania et al. (2010).
	Ceq_CO2 = patch.get_climate().co2 / MM2_PER_M2 / henry_coeff_CO2; // mol L-1
	Ceq_CO2 *= MM2_PER_M2; // mmol m-3

	double henry_coeff_CH4 = henry_k_CH4 * exp(-1.0 * henry_C_CH4 * (1.0 / (surfT + K2degC) - 1.0 / deg25));
	// [L atm mol-1]

	Ceq_CH4 = pp_CH4 / MM2_PER_M2 / henry_coeff_CH4; // mol L-1
	Ceq_CH4 *= MM2_PER_M2; // mmol m-3

	return true;
}


double Soil::diffuse_gas(double Cgas[NLAYERS], double D[NLAYERS], gastype thisgastype, double Ceq, 
					  double kgas, double Dz_m[NLAYERS], double &dailyDiff) {

	// Generic gas diffusion method that works with CO2, CH4 and O2 (as specified with gastype) 
	// See Wania et al. (2010) - Sec 2.5 - for a full description

	// Called like this (e.g. for O2): 
	// diffuse_gas(O2, D_O2, O2gas, Ceq_O2, k_O2, Dz_metre, dailyO2diffusion);
	
	// Atomic mass of this gas [gC/mol]
	double atomic_mass;

	// Concentration of the dissolved gas in question [mmol m-3]
	double C[NLAYERS]; 
	
	if (thisgastype != O2gas)
		atomic_mass = atomiccmass; // i.e. CO2 or CH4 - 12 gC/mol
	else
		atomic_mass = 1.0; // O2 already in mol layer-1
	
	double initialAmount = 0.0;

	for (int ii=IDX; ii<NLAYERS; ii++) {

		initialAmount += Cgas[ii];

		// CO2 & CH4 - g layer-1 to mmol m-3 
		// O2 - mol layer-1 to mmol m-3
		C[ii] = Cgas[ii] / atomic_mass / total_volume_water[ii] * MMOL_PER_MOL; 
	}

	// Set the BC, i.e. equilibrium gas concentrations in the top layer depending on atmospheric concentrations
	// and using Henry's law. See Wania et al. (2010), Eqs. 4-8
	
	// New surface concentration [mmol m-3]
	double Cnew;
	if (thisgastype == O2gas || thisgastype == CH4gas) { // Could also run for CH4 here below.

		if ((Frac_water[IDX] + Frac_water_belowpwp[IDX]) < water_min) { // Could add a snow restriction
			// No diffusion if there's too little liquid water in the top layer
			Cnew = C[IDX]; // Unchanged surface concentration
			dailyDiff = 0.0;
		} 
		else {
			
			// Analytical solution to determine Csurf - see Wania et al. Sec 2.5
			Cnew = Ceq + (C[IDX] - Ceq) * exp(-kgas / volume_liquid_water[IDX]); // mmol m-3
			dailyDiff = (C[IDX] - Cnew) * atomic_mass * volume_liquid_water[IDX] / MMOL_PER_MOL; // mol layer-1 d-1 (O2)
		}

		C[IDX] = Cnew; // mmol m-3 - The new surface concentration
	}

	double surf_conc = C[IDX];
	int layer0 = IDX;

	// Diffusion of gas today?
	if (dailyDiff!=0) {

		// Determine the timestep to use

		double maxD = -0.01;
		// Maximum diffusivity today below the surface layer [m2 d-1]
		for (int ii=IDX+1; ii<NLAYERS; ii++) { 
			// exclude top layer as we've already considered in the top layer above
			if (D[ii] > maxD) {
				maxD = D[ii]; // [m2 d-1]
			}
		}

		double max_timestep = (0.1 * 0.1)/ (2 * maxD); // [day]

		if (max_timestep < Dt_gas && verbosity >= INFO)
			dprintf("%s%8.4f\n","Max gas diffusion timestep exceeded: ",maxD);

		// Assume a tiny diffusivity if there is very little liquid water in the layer
		for (int ii=IDX; ii<NLAYERS; ii++) {
			if ((Frac_water[ii] + Frac_water_belowpwp[ii]) < water_min) {
				D[ii] = 1e-9; // [m2 d-1]
			}
		}

		double C_init[NLAYERS];
		double C_last[NLAYERS];

		// Determine initial and total concentrations
		double totalConc = 0.0; // mmol m-3

		for (int ii=IDX; ii<NLAYERS; ii++) {
			C_init[ii] = C[ii];
			totalConc += C_init[ii];
		}

		double total_diff = 0.0;
		totalConc = 0.0;
		double ratio = 0.0;

		int cncount = 0;

		// Loop until stable concentrations are found, or until 100 iterations have been performed
		int gasdiffix = int(1/Dt_gas);
		
		bool stable = false;

		do {

			totalConc = 0;
			for (int ii=IDX; ii<NLAYERS; ii++) {
				C_init[ii] = C[ii];
				totalConc += C_init[ii];
			}
			
			// Use the CN scheme
			cnstepgas(layer0, D, Dz_m, surf_conc, Dt_gas, C);

			total_diff = 0.0;
			for (int ii=IDX; ii<NLAYERS; ii++) {
				C_last[ii] = C[ii];
				total_diff += fabs(C_last[ii]-C_init[ii]);
			}

			cncount++;

			ratio = total_diff / totalConc;
			
			// stable if the ratio is < 0.01 and we've passed a third of the max iterations 
			if (ratio < 0.01 && cncount > int(gasdiffix/3))
				stable = true;

		} while (cncount < gasdiffix && !stable); // until max 1% variation between iterations, or 100 loops
	}

	// Unit conversions - allocate concentrations to g layer-1 or mole layer-1
	double finalAmount = 0.0;
	for (int ii=IDX; ii<NLAYERS; ii++) {
		Cgas[ii] = C[ii] * atomic_mass * total_volume_water[ii] / MMOL_PER_MOL;
		finalAmount += Cgas[ii];
	}

	// Return the amount of gas that has diffused INTO the soil
	// CO2 & CH4 - g
	// O2 - mol
	double gasIntoSoil = finalAmount - initialAmount;

	if (thisgastype == CH4gas)
		dailyDiff = gasIntoSoil; // gC/m2/d

	return gasIntoSoil;
}


double Soil::calculate_tiller_areas(double r_frac[NLAYERS], double t_area[NLAYERS]) {

	// Calculate the area of tillers
	// See Wania et al. (2010) - Sec 2.6

	double tiller_density = 0.0; // to return [number m-2]

	double graminoid_leafcmass = 0.0; // kgC/m2
	double graminoid_dphen = 0.0;
	double graminoid_anpp_red = 0.0;

	// Loop through this patch object's vegetation and pick out graminoids

	Vegetation& vegetation=patch.vegetation;

	vegetation.firstobj();
	while (vegetation.isobj) {
	
		Individual& indiv=vegetation.getobj();
	
		// Only allow plant-mediated transport from graminoids that are alive
		if (indiv.alive && indiv.pft.has_aerenchyma) {
			// e.g. C3 graminoids on wetlands
			graminoid_leafcmass += indiv.cmass_leaf;
			graminoid_dphen += indiv.phen;
			graminoid_anpp_red += indiv.anpp;
		}
		vegetation.nextobj();
	}

	graminoid_leafcmass *= G_PER_KG; // kgC m-2 to gC m-2

	tiller_density = graminoid_leafcmass * graminoid_dphen / tiller_weight; // m-2

	double tiller_frac = 0.0;

	for (int ii=IDX; ii<NLAYERS; ii++) {

		if (!negligible(tiller_density)) {
			tiller_frac = tiller_density * r_frac[ii];

			// Consider this as the fraction of the area (1m2) occupied by tillers in this later
			t_area[ii] = PI * pow(tiller_radius, 2.0) * tiller_frac * tiller_por; // [unitless m2 m-2]

			if (t_area[ii] < 0.000000001)
				t_area[ii] = 0.0;
		}
		else {
			t_area[ii] = 0.0;
		}
	}

	return tiller_density;
}


bool Soil::plant_gas_transport(double Cgas[NLAYERS], double Ceq, double kgas, gastype thisgastype, 
							 double& plantTransportToday) {

	// Plant transport of O2 or CH4
	// See Wania et al. (2010) - Sec 2.6

	// The limiting factor to plant transport are the number of tillers as
    // they represent the cross-sectional area available to gas transport.
	// The biomass is related to tiller density based on Schimel 1995. The
	// average biomass was 185g/m2 and the average number of tillers was 380/m2.
    // This gives us 2.05 tillers per g biomass and 0.48 g per tiller, which
    // corresponds to 0.22 g C per tiller


	double Cnew = 0.0; // mmol m-3
	double atomic_mass;

	if (thisgastype != O2gas)
		atomic_mass = atomiccmass; // i.e. CO2 or CH4 - 12 g/mol
	else
		atomic_mass = 1.0; 

	// to return
	plantTransportToday = 0.0;

	for (int ii=IDX; ii<NLAYERS; ii++) {

		// Weight plant transport by the area of porous root cross-sections. 
		// plant_trans(i) [mmol layer-1 d-1]
        // the 2.07 are from the equation for k_CH4, when U10 is zero.

		Cgas[ii] = Cgas[ii] / atomic_mass / total_volume_water[ii] * MMOL_PER_MOL;
		// CH4 and CO2: gC layer-1 to mmol m-3
		// O2: mol layer-1 to mmol m-3

		if ((Frac_water[ii] + Frac_water_belowpwp[ii]) < water_min || negligible(tiller_area[ii]) || volume_liquid_water[ii] < 0.0001) {
			Cnew = Cgas[ii]; // no plant transport to/from this layer
		} 
		else {
			double exponent = kgas / volume_liquid_water[ii];
			exponent *= tiller_area[ii];
			double fac = exp(-1.0 * exponent);
			double change = Cnew + (Cgas[ii] - Ceq) * fac;

			Cnew = Ceq + (Cgas[ii] - Ceq) * exp(-kgas / (volume_liquid_water[ii] / tiller_area[ii])); // mmol m-3
			plantTransportToday += (Cgas[ii] - Cnew) 
				* atomic_mass * total_volume_water[ii] / MMOL_PER_MOL; // mol layer-1 d-1 (O2) or // gC m-2 d-1 (CH4 & CO2)

			// Restrict TINY negative plant transport
			if (plantTransportToday < 0.0 && fabs(plantTransportToday) < 1e-8) {
				plantTransportToday = 0.0;
				Cnew = Cgas[ii];
			}
		}
		
		// mmol m-3 to gC layer-1 (CO2 or CH4), or mol layer-1 (O2) 
		Cgas[ii] = Cnew * atomic_mass * total_volume_water[ii] / MMOL_PER_MOL;
	}
	return true;
}


bool Soil::calculate_gas_ebullition(double& ebull_today) {

	// Gas ebullition
	// Returns the gas ebullition today
	// See Wania et al. (2010) - Sec 2.7

	double waterheight = 0.0;
	ebull_today = 0.0;

	double wtpp = wtp[date.day];	// Today's water table position (-300 to +100) mm
	bool emitToAtmosphere = true;	// Whether to emit CH4 directly to the atmosphere (true by default)
	int bubbleToLayer = IDX;		// The soil layer to bubble to.

	if (wtpp <= 0) {
		
		if (wtpp >= -100.0)
			bubbleToLayer = IDX;
		else if (wtpp >= -200.0 && wtpp < -100.0)
			bubbleToLayer = IDX+1;
		else 
			bubbleToLayer = IDX+2;
	} 
	else { // Standing water
		bubbleToLayer = IDX-1;
	}

	// Ensure ebullition from all layers
	bubbleToLayer = IDX-1;

	for (int ii=IDX; ii<NLAYERS; ii++) {

		double TsoilK = T_soil[ii] + K2degC;

		waterheight += total_volume_water[ii]; // m

		// Bubble formation below the WT only
		if (ii > bubbleToLayer) {

			// Max CH4 that can be dissolved (Wania et al. (2010), Eqn 15)
			double CH4_diss_max = 0.05708 - 0.001545 * max(0.0, T_soil[ii]) + 
				0.00002069 * max(0.0, pow(T_soil[ii],2.0)); // ml CH4 ml-1 H2O

			// Water pressure
			double hydro_press = rho_H2O * gravity * waterheight; // N, or kg m s-2

			// Maximum volume that can be dissolved
			double CH4_diss_max_m3 = CH4_diss_max * volume_liquid_water[ii]; // m3 CH4 layer-1

			// Use ideal gas law to convert to mol CH4 layer-1
			double CH4_diss_max_mol = CH4_diss_max_m3 * (hydro_press + atm_press) / (R_gas * TsoilK); // mol CH4 layer-1

			// Convert to g CH4-C layer-1
			double CH4_diss_max_g = CH4_diss_max_mol * mr_C; // g CH4-C layer-1
			
			CH4_ebull_ind[ii] = 0.0;
			
			// Restrict ebullition to cases when there is enough liquid water and when soil T > 0. 
			if ((Frac_water[ii] + Frac_water_belowpwp[ii]) > water_min && T_soil[ii] > 0.0) {

				double henry_k_cc_CH4 = TsoilK / (12.2 * henry_k_CH4);
				CH4_diss[ii] = min(CH4_diss_max_g, henry_k_cc_CH4 * CH4[ii]);
				
				// Override (a la Wania et al. 2010) with this new, simpler ebullition
				CH4_diss[ii] = min(CH4_diss_max_g, CH4[ii]);

				CH4_gas[ii] = CH4[ii] - CH4_diss[ii]; // g CH4-C layer-1 - gas
				CH4_gas_vol[ii] = CH4_gas[ii] / atomiccmass * R_gas * TsoilK / (hydro_press + atm_press); // [m3]
				CH4_vgc[ii] = CH4_gas_vol[ii] / (Dz[ii] / MM_PER_M * SQ_M); // m3/m3
			
				// Ebullition/bubble formation if the volumetric gas content (VGC) exceeds vgc_high * bubble_CH4_frac, 
				// where bubble_CH4_frac id the CH4 fraction of gas bubbles
				if (CH4_vgc[ii] > vgc_high * bubble_CH4_frac) {
				
					CH4_ebull_vol[ii] = (CH4_vgc[ii] - vgc_low * bubble_CH4_frac) * Dz[ii] / MM_PER_M * SQ_M; // m3/m3
					CH4_vgc[ii] = vgc_low * bubble_CH4_frac; 
					CH4_ebull_ind[ii] = CH4_ebull_vol[ii] * (hydro_press + atm_press) / (R_gas * TsoilK); // mol
					CH4_gas[ii]	= CH4_vgc[ii] * Dz[ii] / MM_PER_M * SQ_M * (hydro_press + atm_press) / (R_gas * TsoilK); // mol
					CH4_ebull_ind[ii] *= atomiccmass; // mol to g
					CH4_gas[ii] *= atomiccmass; // mol to g 

					// Update the amount of CH4 to emit today, from this layer
					ebull_today += CH4_ebull_ind[ii]; // g m-2 d-1
				}

				// Now recalculate the total CH4 amount in this layer, both dissolved and gaseous
				CH4[ii] = CH4_diss[ii] + CH4_gas[ii];

			} 
			else {
			
				// No change in the dissolved CH4 amount
				CH4_diss[ii] = CH4_diss_yesterday[ii]; 
			} // Frac_water
		} 
		else {
			
			// No change in the dissolved CH4 amount
			CH4_diss[ii] = CH4_diss_yesterday[ii]; 	
		} // bubble layer check
	} // for ii

	// If the water table is below the surface, then put the bubbled methane into the first unsaturated layer
	if (!emitToAtmosphere) {

		CH4[bubbleToLayer] += ebull_today;  
		ebull_today = 0.0; // No emission today.

		// O2 is in [mol layer-1]
		O2[bubbleToLayer] *= oxid_frac; // since 25% of O2 used by roots themselves...  
		CH4[bubbleToLayer] /= atomiccmass; // mol layer-1
		CH4_oxid[bubbleToLayer] = min(CH4[bubbleToLayer], 0.5 * O2[bubbleToLayer]); // usually 75%
		CH4[bubbleToLayer] = (CH4[bubbleToLayer] - CH4_oxid[bubbleToLayer]) * atomiccmass; // gC layer again	
		O2[bubbleToLayer] -= 2.0 * CH4_oxid[bubbleToLayer]; // subtract the moles used in oxidation

		CO2_soil[bubbleToLayer] += CH4_oxid[bubbleToLayer] * atomiccmass; // Oxidised CH4 becomes CO2 (and water!)
	}

	return true;
}


double Soil::get_ch4_content() {
	// Return dissolved CH4 content. Units: g CH4-C / m2
	return ch4_store;
}

double Soil::get_co2_content() {
	// Return dissolved CO2 content. Units: g C / m2
	return co2_store;
}

void Soil::calculate_carbon_store(int dy, bool today) {

	// Carbon accounting routine - updates co2_store and ch4_store 
	// See Wania et al. (2010)

	int day = 0;

	// Reset the stores
	co2_store = 0.0;
	ch4_store = 0.0;

	if (today)
		day = dy;
	else
		day = dy-1;

	if (dy > 0 || today) {
		for (int ii=IDX; ii<NLAYERS; ii++) {
			co2_store += CO2_soil[ii];
			ch4_store += CH4[ii];
		}
	} 
	else {
		// Jan 1 or yesterday
		for (int ii=IDX; ii<NLAYERS; ii++) {
			co2_store += CO2_soil_yesterday[ii];
			ch4_store += CH4_yesterday[ii];
		}
	}
}


bool Soil::methane(bool generatemethane) {

	// Main methane routine, called daily. 
	// Calculates daily methane fluxes and updates methane concentrations in peatland soil layers.
	// Algorithm etc. from Wania et al. 2010
	// This version coded by Paul Miller, based on Rita Wania's F90 code

	if (!generatemethane) {

		// No generation of methane for what ever reason, so report the heterotrophic respiration 
		// flux calculated from the call to som_dynamics(), and 0 for the 4 CH4 fluxes.
		patch.fluxes.report_flux(Fluxes::SOILC, dcflux_soil); // kgC/m2
		patch.fluxes.report_flux(Fluxes::CH4C, 0.0);
		patch.fluxes.report_flux(Fluxes::CH4C_DIFF, 0.0);
		patch.fluxes.report_flux(Fluxes::CH4C_PLAN, 0.0);
		patch.fluxes.report_flux(Fluxes::CH4C_EBUL, 0.0);

	} 
	else if (patch.get_climate().lat < PEATLAND_WETLAND_LATITUDE_LIMIT) {
		
		// Simple Spahni et al. (2011) approach in which a fraction (CH4toCO2_inundated) of the daily 
		// heterotrophic respiration is assumed to be in the form of CH4

		double inundated_CH4_flux_today = dcflux_soil * G_PER_KG * CH4toCO2_inundated; // g CH4-C/m2/day

		// Report the heterotrophic respiration and the 4 CH4 fluxes.
		patch.fluxes.report_flux(Fluxes::SOILC, dcflux_soil - inundated_CH4_flux_today * KG_PER_G); // kgC/m2
		patch.fluxes.report_flux(Fluxes::CH4C, inundated_CH4_flux_today); // g CH4-C/m2/day
		patch.fluxes.report_flux(Fluxes::CH4C_DIFF, 0.0);
		patch.fluxes.report_flux(Fluxes::CH4C_PLAN, 0.0);
		patch.fluxes.report_flux(Fluxes::CH4C_EBUL, 0.0);

	} 
	else {
        // Max allowed error in checks
        const double MAX_ERR = 0.000001;
		const double MAX_ERR_BALANCE = 0.0001;
		const double LARGE_ERR = 0.01;

		// variables for debugging
		bool allow_planttransport = true;
		bool allow_ebullition = true;
		bool allow_ch4diffusion = true;
		bool allow_o2diffusion = true;

		// Layer diffusivities
		double D_CH4[NLAYERS];
		double D_CO2[NLAYERS];
		double D_O2[NLAYERS];

		double Dz_metre[NLAYERS]; // layer depths [m]

		// Daily diffusion 
		double O2_diff_today;
		double CO2_diff_today;
		double CH4_diff_today;
	
		// Daily plant transport
		double O2_plant_today;
		double CH4_plant_today;
		double CO2_plant_today;
	
		// Total daily ebullition [g CH4-C d-1]
		double CH4_ebull_today;
		
		int daynum = date.day;

		// Initialise the peatland root fractions.
		if (daynum == 0)
			init_peatland_root_fractions();

		// Budget/conservation variables [gC/m2]
		double ch4_c_in = 0.0;
		double co2_c_in = 0.0;

		double ch4_c_store_init = 0.0;
		double co2_c_store_init = 0.0;

		double ch4_c_store_now = 0.0;
		double co2_c_store_now = 0.0;

		// Initialise components of total_C_flux, CO2_flux and CH4_flux
		O2_diff_today = 0.0;
		CO2_diff_today = 0.0;
		CH4_diff_today = 0.0;
		CH4_ebull_today = 0.0;
		CH4_plant_today = 0.0;
		CO2_plant_today = 0.0;
		O2_plant_today = 0.0;

		// *** STEP 0 ***

		// C budget before today's methane calculations

		calculate_carbon_store(daynum, false);
		ch4_c_store_init = ch4_store;
		co2_c_store_init = co2_store;

		// *** STEP 1 ***
		// Diffusivities and layer depths, and update of gas constants

		// Set layer gas diffusivities 
		calculate_gas_diffusivities(D_CH4, D_CO2, D_O2);
		// Units: m2 d-1

		// Layer depths in m, and volume of water in m3

		for (int ii=IDX; ii<NLAYERS; ii++) {
			Dz_metre[ii] = Dz[ii] / MM_PER_M; // layer depths in m
			volume_liquid_water[ii] = (Frac_water[ii] + Frac_water_belowpwp[ii]) * Dz_metre[ii];	// m3 water
			total_volume_water[ii] = (Frac_water[ii] + Fpwp_ref[ii] + Frac_ice[ii]) * Dz_metre[ii];	// m3 water + ice
		}

		// If we have standing water...
		if (wtp[daynum] > 0.0) {
			// wtp in mm
			Dz_metre[IDX] += wtp[daynum] / MM_PER_M;
			volume_liquid_water[IDX] += (Frac_water[MIDX] + Frac_water_belowpwp[MIDX]) * wtp[daynum] / MM_PER_M;
			total_volume_water[IDX] += (Frac_water[MIDX] + Fpwp_ref[MIDX] + Frac_ice[MIDX]) * wtp[daynum] / MM_PER_M;
		}

		// Update the temperature-dependent gas parameters
		bool gasParamsOK = update_daily_gas_parameters();
		// Units: 
		// k_O2 etc: [m d-1]
		// Ceq_O2: [mmol m-3]

		if (!gasParamsOK)
			return false;

		// *** STEP 2 ***
		// Calculate CH4 & CO2 production	

		// Daily decomposition, converted to gC/m2 from kgC/m2
		double drh = dcflux_soil * G_PER_KG;
		// Today's decomposition of litter and soil C pools

		double c_input = 0.0; // conservation check - C in methane produced
		double ch4_init_prod = 0.0;
		double total_ch4 = 0.0;

		for (int ii=IDX; ii<NLAYERS; ii++) {

			double anoxic = 1.0 - Frac_air[ii] - Fgas;
			double c_input_old = c_input;

			// *** CH4 PRODUCTION IN THIS LAYER ***

			if ((Frac_water[ii] + Frac_water_belowpwp[ii]) < water_min)
				CH4_prod[ii] = 0.0;
			else
				CH4_prod[ii] = anoxic * CH4toCO2_peat * rootfrac[ii] * drh; // gC/m2
	
			c_input += CH4_prod[ii];
			ch4_init_prod += CH4_prod[ii];

			// Add CH4 production to CH4 pool 
			CH4[ii] = CH4_yesterday[ii] + CH4_prod[ii];
			total_ch4 += CH4[ii];

			// *** CO2 PRODUCTION ***

			CO2_soil_prod[ii] = rootfrac[ii] * drh - CH4_prod[ii];
			c_input += CO2_soil_prod[ii];

			// CO2 pool
			// Add CO2 production to CO2 pool 
			CO2_soil[ii] = CO2_soil_yesterday[ii] + CO2_soil_prod[ii];

			double c_added = c_input - c_input_old;
			double ought_to_have_been = rootfrac[ii] * drh;

			double er = ought_to_have_been - c_added;
		}

		// C conservation test:
		if (fabs(drh - c_input) > MAX_ERR && verbosity >= WARNING) {
			dprintf("%s%8.5f\n","Bad C conservation at the start of Soil::methane()",fabs(drh - c_input));	
			//return false;
		}

		// debugging:
		if ((total_ch4 < 0.0000000 || total_ch4 > 10000000) && verbosity >= WARNING) {
			dprintf("%s%8.5f\n","Bad total_ch4 at the start of Soil::methane()",total_ch4);	
		}

		// C conservation test:
		calculate_carbon_store(daynum, true);
		ch4_c_store_now = ch4_store;
		co2_c_store_now = co2_store;
		double check = ch4_c_store_now + co2_c_store_now - ch4_c_store_init - co2_c_store_init - drh;
		if (fabs(check) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","C imbalance in Soil::methane() : After C in: ",check);

		// *** STEP 3 ***
		// Diffusion of O2

		double dailyO2diffusion = 0.0;
		double molesO2IntoSoil = 0.0; 
		if (allow_o2diffusion && dsnowdepth < 50.0) // no O2 diffusion until snow depth < 50mm 
			molesO2IntoSoil = diffuse_gas(O2, D_O2, O2gas, Ceq_O2, k_O2, Dz_metre, dailyO2diffusion);
		O2_diff_today = dailyO2diffusion; // mol O2 into the soil (and then diffused downwards) 
		// Should be negative, i.e. O2 diffuses INTO the soil

		// *** STEP 4 ***
		// Plant transport of oxygen

		// Tiller set-up
		// Loop through the individuals present and sum the leaf carbon mass for C3-graminoids 
		double tillers = calculate_tiller_areas(rootfrac, tiller_area); // m-2
	
		// O2 plant transport
		double plantTransportToday_O2 = 0.0;
	
		if (allow_planttransport && tillers > 0) // Only allow plant transport if there are graminoid tillers present
			plant_gas_transport(O2, Ceq_O2, k_O2, O2gas, plantTransportToday_O2);

		O2_plant_today = plantTransportToday_O2; // mol O2 into soil through plants
	
		// *** STEP 5 ***
		// Diffusion of CH4

		double ch4_before_diffusion = ch4_store;

		double dailyCH4diffusion = 0.0;
		double gramCH4IntoSoil = 0.0; 
		if (allow_ch4diffusion && dsnowdepth < 50.0) // no CH4 diffusion until snow depth < 50mm 
			gramCH4IntoSoil = diffuse_gas(CH4, D_CH4, CH4gas, Ceq_CH4, k_CH4, Dz_metre, dailyCH4diffusion);
		CH4_diff_today = -dailyCH4diffusion; // SHOULD BE 0 - gC m-2 d-1
		// Should be positive, i.e. upward flux

		if (fabs(CH4_diff_today) < MAX_ERR)
			CH4_diff_today = 0.0; // remove tiny values

		// C conservation test:
		if ((CH4_diff_today < -0.1 || CH4_diff_today > 10000000) && verbosity >= WARNING) {
			dprintf("%s%8.5f\n","Bad CH4 diffusion in Soil::methane()",CH4_diff_today);	
		}

		// C conserved?
		calculate_carbon_store(daynum, true);
		check = ch4_store -ch4_before_diffusion+CH4_diff_today;
		if (fabs(check) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after diffusion of CH4: ",check);

		// *** STEP 6 ***
		// CH4 oxidation

		// CH4 units: gC layer-1 
		// O2 units: mol layer-1

		// Assume that 1/2 of the O2 is utilized by other electron
		// acceptors and you need 2 moles of O2 to oxidise 1 mole of CH4.
		for (int ii=IDX; ii<NLAYERS; ii++) {
			// O2 is in [mol layer-1]
			O2[ii] *= oxid_frac; // since ((1-oxid_frac)*100)% of O2 used by roots themselves... [mol layer-1] 
			CH4[ii] /= atomiccmass; // mol layer-1
			CH4_oxid[ii] = min(CH4[ii], 0.5 * O2[ii]); // usually 75%
			CH4[ii] = (CH4[ii] - CH4_oxid[ii]) * atomiccmass; // gC layer again	
			O2[ii] -= 2.0 * CH4_oxid[ii]; // subtract the moles used in oxidation

			CO2_soil[ii] += CH4_oxid[ii] * atomiccmass; // Oxidised CH4 becomes CO2 (and water!)
		}

		// *** STEP 6b ***
		// CH4 diffusion from top layer AFTER diffusion and oxidation
		// Now done after CN diffusion

		double ch4_diff_before_BC = CH4_diff_today;

		// Check C budget
		double ch4_before_oxid_and_diffusion = ch4_store;
		double co2_before_oxid_and_diffusion = co2_store;
	
		double checkbefore = ch4_before_oxid_and_diffusion + co2_before_oxid_and_diffusion + ch4_diff_before_BC;

		calculate_carbon_store(daynum, true);
		check = ch4_store+co2_store+CH4_diff_today;

		double checkafter1 = check - checkbefore;

		if (fabs(checkafter1) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after oxidation and diffusion of CH4: ",checkafter1);

		checkbefore = check;

		// *** STEP 7 ***
		// Plant transport of CH4

		double plantTransportToday_CH4 = 0.0;

		if (allow_planttransport && tillers > 0) // Only allow plant transport of there are gramoinoid tillers present
			plant_gas_transport(CH4, Ceq_CH4, k_CH4, CH4gas, plantTransportToday_CH4);

		CH4_plant_today = plantTransportToday_CH4; // CH4 into soil through plants (gC m-2)

		// Check C budget
		calculate_carbon_store(daynum, true);
		check = ch4_store+co2_store+CH4_diff_today+CH4_plant_today;
		double checkafter2 = check - checkbefore;
		if (fabs(checkafter2) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after plant transport of CH4: ",checkafter2);
		checkbefore = check;

		// *** STEP 8 ***
		// Ebullition of CH4
		if (allow_ebullition)
			calculate_gas_ebullition(CH4_ebull_today);
		// reduces CH4 and updates CH4_ebull_today

		// Check C budget
		calculate_carbon_store(daynum, true);
		check = ch4_store+co2_store+CH4_diff_today+CH4_plant_today+CH4_ebull_today;

		double checkafter3 = check - checkbefore;

		if (fabs(checkafter3) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after ebullition of CH4: ",checkafter3);

		// STEPS 9 & 10
		// Instead of calling diffusion and plant transport for CO2, I simply diffuse ALL the CO2. 
		// This should improve C conservation
		// In future one could look at how much IS emitted, especially as CO2 could be saved and emitted in 
		// bursts during spring thaw

		CO2_diff_today = co2_store; 

		// Correct & record Dec 31 data
		for (int ii=IDX; ii<NLAYERS; ii++) {
			CO2_soil[ii] = 0.0;
		}
		
		// Recalculate ch4_store and co2_store now, after CO2 diffusion
		calculate_carbon_store(daynum, true);

		CO2_plant_today = 0.0;

		// *** STEP 9 ***
		// Conservation steps

		// Total C flux out of the soil today [gC m-2 d-1] 
		double total_C_flux = CO2_diff_today + CH4_diff_today + CH4_ebull_today + CH4_plant_today + CO2_plant_today;

		// CO2_flux_today = CO2_diff_today + CO2_plant_today;
		// Calculate daily CH4 flux [gCH4-C m-2 d-1]
		double CH4_flux_today = CH4_diff_today + CH4_plant_today + CH4_ebull_today;

		// Reduce heterotrophic respiration by this CH4-C amount

		// Report the heterotrophic respiration and the 4 CH4 fluxes.
		patch.fluxes.report_flux(Fluxes::SOILC, dcflux_soil-CH4_flux_today/G_PER_KG); // Units for Fluxes::SOILC are kgC/m2/time
		patch.fluxes.report_flux(Fluxes::CH4C, CH4_flux_today); // gC/m2
		patch.fluxes.report_flux(Fluxes::CH4C_DIFF, CH4_diff_today);
		patch.fluxes.report_flux(Fluxes::CH4C_PLAN, CH4_plant_today);
		patch.fluxes.report_flux(Fluxes::CH4C_EBUL, CH4_ebull_today);

		// Total C content of soil
		double c_soil = 0.0;
		double co2_soil = 0.0;

		calculate_carbon_store(daynum, true);
		c_soil += ch4_store + co2_store;

		if (total_C_flux < -LARGE_ERR && verbosity >= WARNING) {
			dprintf("%s%10.5f\n","Large negative C flux (total_C_flux) in Soil::methane()",total_C_flux);	
		}

		double final_C_budget = (ch4_c_store_init + co2_c_store_init + c_input) - (c_soil + total_C_flux);
		if (fabs(final_C_budget) > MAX_ERR_BALANCE && verbosity >= WARNING) {
			dprintf("%s%10.5f\n","Unbalanced C budget in Soil::methane()",total_C_flux);	
		}

		// Record today's values
		for (int ii=IDX; ii<NLAYERS; ii++) {
			CO2_soil_yesterday[ii] = CO2_soil[ii]; 
			CH4_yesterday[ii] = CH4[ii]; 
			CH4_diss_yesterday[ii] = CH4_diss[ii];
			CH4_gas_yesterday[ii] = CH4_gas[ii];
		}
	} // generatemethane?

	return true;
}



void Soil::init_peatland_root_fractions() {

	// Initialise the root fractions in each layer of the wetland, assuming an exponential decrease of root biomass with depth
	// in accordance with observations from peat cores.

	// Called on the first day of the simulaton.
	
	// See Wania et al. (2010) - Section 2.2.2

	const double b = 25.1695; // cm
	const double corr = 2.49994026;

	double layermiddepth = Dz_soil/(2 * MM_PER_CM); // normally 5 cm

	double sumrootfrac = 0.0;

	for (int ii=IDX; ii<NLAYERS; ii++) {

		layermiddepth -= 10.0; 
		
		if (ii == NLAYERS-1)
			rootfrac[ii] = 1.0 - sumrootfrac;
		else
			rootfrac[ii] = exp(layermiddepth/b) / corr;		

		sumrootfrac += rootfrac[ii];
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Wania, R., Ross, I., and Prentice, I. C.: Integrating peatlands and permafrost into a dynamic
//   global vegetation model; 1, Evaluation and sensitivity of physical land surface processes, Global
//   Biogeochemical Cycles, 23, doi:10.1029/2008gb003412, 2009a.
// Wania, R., Ross, I., and Prentice, I. C.: Integrating peatlands and permafrost into a dynamic
//   global vegetation model; 2, Evaluation and sensitivity of vegetation and carbon cycle processes,
//   Global Biogeochemical Cycles, 23, doi:10.1029/2008gb003413, 2009b.
// Wania, R., Ross, I., and Prentice, I. C.: Implementation and evaluation of a new methane model
//   within a dynamic global vegetation model: LPJ-WHyMe v1.3.1, Geosci Model Dev, 3, 565-584,
//   doi:DOI 10.5194/gmd-3-565-2010, 2010. 
