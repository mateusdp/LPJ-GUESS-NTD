///////////////////////////////////////////////////////////////////////////////////////
/// \file emdi.h
/// \brief Extra code used by the EMDI benchmarks
///
/// \author Joe Siltberg
/// $Date: 2013-10-10 10:20:33 +0200 (Do, 10. Okt 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_EMDI_H
#define LPJ_GUESS_EMDI_H

#include <gutil.h>
#include <map>
#include <utility>
#include "cruinput.h"
#include "guess.h"

/// The EMDI input module
/** This is a subclass of the CRU input module, and works
 *  exactly the same except that plant available water content
 *  is overridden for each gridcell with values taken from
 *  the gridlist.
 */
class EMDIInput : public CRUInput {
public:

	void init();

	bool getgridcell(Gridcell& gridcell);

private:

	void rememberPAWC(double dlon, double dlat, xtring pawc);
	
	void overrideAWC(double lon, double lat, Soiltype& soiltype);

	std::map<std::pair<double, double>, double> pawcPerGridCell;
};

#endif // LPJ_GUESS_EMDI_H
