///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete.h
/// \brief Distribution of N among individuals according to supply, demand and
///        the individuals' uptake strength
///
/// \author David Wårlind
/// $Date: 2013-04-30 11:05:01 +0200 (Tue, 30 Apr 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_NCOMPETE_H
#define LPJ_GUESS_NCOMPETE_H

#include <vector>

/// Represents an individual competing for nitrogen uptake
/** Contains what the ncompete function below needs to know about
 *  an individual in order to do the distribution of N.
 */
struct NCompetingIndividual {
	/// This individual's N demand
	double ndemand;

	/// A meassure of this individual's uptake strength
	double strength;

	/// Output from ncompete - fraction of the demand satisfied by the distribution
	double fnuptake;
};

/// Distributes N among individuals according to supply, demand and uptake strength
/** Grasses should get at least 5% and no
 *  individual should get more than 100% of its nitrogen demand. 
 */
void ncompete(std::vector<NCompetingIndividual>& individuals, 
              double nmass_avail);

/// Represents an individual competing for phosphorus uptake
/** Contains what the pcompete function below needs to know about
*  an individual in order to do the distribution of P.
*/
struct PCompetingIndividual {
	/// This individual's P demand
	double pdemand;

	/// A meassure of this individual's uptake strength
	double strength;

	/// Output from ncompete - fraction of the demand satisfied by the distribution
	double fpuptake;
};

/// Distributes P among individuals according to supply, demand and uptake strength
/** Grasses should get at least 5% and no
*  individual should get more than 100% of its phosphorus demand.
*/
void pcompete(std::vector<PCompetingIndividual>& individuals,
	double pmass_avail);

#endif // LPJ_GUESS_NCOMPETE_H
