///////////////////////////////////////////////////////////////////////////////////////
/// \file cru.cpp
/// \brief Functions for reading the CRU-NCEP data set
///
/// $Date: 2013-11-04 16:30:55 +0100 (Mon, 04 Nov 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "cru_ts30.h"
#include <stdio.h>
#include <math.h>
#include <vector>

// header files for the CRU-NCEP data archives
#include "cruncep_1901_2015.h"
#include "cruncep_1901_2015misc.h"

namespace CRU_TS30 {
 
bool searchcru(char* cruark,double dlon,double dlat,int& soilcode,
               double mtemp[NYEAR_HIST][12],
               double mprec[NYEAR_HIST][12],
               double msun[NYEAR_HIST][12]) {

	// !!!! NEW VERSION OF THIS FUNCTION - guess2008 - NEW VERSION OF THIS FUNCTION !!!!
	// Please note the new function signature. 

	// Archive object. Definition in new header file, cru.h
	Cruncep_1901_2015Archive ark;

	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cruncep_1901_2015 data; // struct to hold the data

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon;
		data.lat = dlat;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays. 
		soilcode=(int)data.soilcode[0];


		for (y = 0; y < NYEAR_HIST; y++) {
			for (m=0;m<12;m++) {
				mtemp[y][m] = data.mtemp[y*12+m]; // now degC
				mprec[y][m] = data.mprec[y*12+m]; // mm (sum over month)
				
				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (mprec[y][m] <= 1.0) mprec[y][m] = 0.0;
				
				msun[y][m]  = data.mswrad[y*12+m];   // shortwave radiation

			}
		}


		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}



bool searchcru_misc(char* cruark,double dlon,double dlat,int& elevation,
                    double mfrs[NYEAR_HIST][12],
                    double mwet[NYEAR_HIST][12],
                    double mdtr[NYEAR_HIST][12]) {
	
	// Please note the new function signature. 

	// Archive object
	Cruncep_1901_2015miscArchive ark; 
	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cruncep_1901_2015misc data;

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon;
		data.lat = dlat;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays.
		// Note that the multipliers are NOT the same as in searchcru above!
		elevation=(int)data.elv[0]; // km * 1000?

		for (y = 0; y < NYEAR_HIST; y++) { 
			for (m=0;m<12;m++) {

				// guess2008 - catch rounding errors 
				mfrs[y][m] = data.mfrs[y*12+m]; // days
				if (mfrs[y][m] < 0.1) 
					mfrs[y][m] = 0.0; // Catches rounding errors

				mwet[y][m] = data.mwet[y*12+m]; // days
				if (mwet[y][m] <= 0.1) 
					mwet[y][m] = 0.0; // Catches rounding errors

				mdtr[y][m] = data.mdtr[y*12+m];  // degC

				// For some reason there are negative dtr values in
				// the CRU binaries(!). Set these to zero for now.
				mdtr[y][m] = max(0.0, mdtr[y][m]);

				/*
				If vapour pressure is needed:
				mvap[y][m] = data.mvap[y*12+m];
				*/
			}
		}

		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}


bool findnearestCRUdata(double searchradius, char* cruark, double& lon, double& lat, 
                        int& scode, 
                        double hist_mtemp1[NYEAR_HIST][12], 
                        double hist_mprec1[NYEAR_HIST][12], 
                        double hist_msun1[NYEAR_HIST][12]) {

	// First try the exact coordinate
	if (searchcru(cruark, lon, lat, scode, hist_mtemp1, hist_mprec1, hist_msun1)) {
		return true;
	}
	
	if (searchradius == 0) {
		// Don't try to search
		return false;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(lon*2)/2 + 0.25;
	double center_lat = floor(lat*2)/2 + 0.25;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy 
	// sorting.
	using std::pair;
	using std::make_pair;
	typedef pair<double, double> point;
	std::vector<pair<double, point> > search_points;

	const double STEP = 0.5;
	const double EPS = 1e-15;

	for (double y = center_lon-searchradius; y <= center_lon+searchradius+EPS; y += STEP) {
		for (double x = center_lat-searchradius; x <= center_lat+searchradius+EPS; x += STEP) {
			double xdist = x - lat;
			double ydist = y - lon;
			double dist = sqrt(xdist*xdist + ydist*ydist);
			
			if (dist <= searchradius + EPS) {
				search_points.push_back(make_pair(dist, make_pair(y, x)));
			}
		}
	}

	// Sort by increasing distance
	std::sort(search_points.begin(), search_points.end());

	// Find closest coordinate which can be found in CRU
	for (unsigned int i = 0; i < search_points.size(); i++) {
		point search_point = search_points[i].second;
		double search_lon = search_point.first;
		double search_lat = search_point.second;

		if (searchcru(cruark, search_lon, search_lat, scode, 
		              hist_mtemp1, hist_mprec1, hist_msun1)) {
			lon = search_lon;
			lat = search_lat;
			return true;
		}
	}

	return false;
}

}
