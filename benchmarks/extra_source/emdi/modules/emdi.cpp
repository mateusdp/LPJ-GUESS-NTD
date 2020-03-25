///////////////////////////////////////////////////////////////////////////////////////
/// \file emdi.cpp
/// \brief Extra code used by the EMDI benchmarks
///
/// \author Joe Siltberg
/// $Date: 2013-07-04 14:55:26 +0200 (Do, 04. Jul 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "emdi.h"
#include <sstream>

REGISTER_INPUT_MODULE("emdi", EMDIInput);

void EMDIInput::init() {
	CRUInput::init();

	// Read in the gridlist again to get plant available water content

	bool eof = false;

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	while (!eof) {

		double dlon, dlat;
		xtring descrip;

		// Read next record in file
		eof=!readfor(in_grid,"f,f,a#",&dlon,&dlat,&descrip);

		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			rememberPAWC(dlon, dlat, descrip);
		}
	}


	fclose(in_grid);
	
}

bool EMDIInput::getgridcell(Gridcell& gridcell) {
	if (CRUInput::getgridcell(gridcell)) {
		overrideAWC(gridcell.get_lon(), gridcell.get_lat(), gridcell.soiltype);
		return true;
	}
	else {
		return false;
	}
}

void EMDIInput::rememberPAWC(double dlon, double dlat, xtring pawc) {
	std::istringstream is((char*)pawc);
	double first_number_in_string;
	is >> first_number_in_string;
		
	pawcPerGridCell[std::make_pair(dlon, dlat)] = first_number_in_string;
}
	
void EMDIInput::overrideAWC(double lon, double lat, Soiltype& soiltype) {
	double pawc = pawcPerGridCell[std::make_pair(lon, lat)];
	if (pawc > 0) {
		// Assume pawc applies to the whole 1.5m, so replace soiltype.awc with a scaled pawc
		double scaled = pawc / 300.0; // as pawc applies to the upper 30cm
		soiltype.awc[0]=SOILDEPTH_UPPER*scaled;
		soiltype.awc[1]=SOILDEPTH_LOWER*scaled;
	}
}
