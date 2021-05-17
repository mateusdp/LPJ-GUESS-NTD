///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete_test.cpp
/// \brief Unit tests for nitrogen uptake competition
///
/// \author Joe Siltberg
/// $Date: 2014-06-16 15:46:53 +0200 (Mon, 16 Jun 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "catch.hpp"

#include "ncompete.h"

TEST_CASE("ncompete/single", "Testing a single individual") {
	std::vector<NCompetingIndividual> indivs(1);

	indivs[0].ndemand = 1;
	indivs[0].strength = 1;

	// more N available than needed
	ncompete(indivs, 2.0);

	REQUIRE(indivs[0].fnuptake == Approx(1));


	// less N available than needed
	ncompete(indivs, 0.5);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));

	// more N available than needed
	ncompete(indivs, 2.0);

	REQUIRE(indivs[0].fnuptake == Approx(1));


	// less N available than needed
	ncompete(indivs, 0.5);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));	
}


TEST_CASE("ncompete/double", "Testing two individuals") {
	std::vector<NCompetingIndividual> indivs(2);

	// two equal individuals

	indivs[0].ndemand = 1;
	indivs[0].strength = 1;
	
	indivs[1] = indivs[0];

	// more N available than needed
	ncompete(indivs, 3.0);
	
	REQUIRE(indivs[0].fnuptake == Approx(1));
	REQUIRE(indivs[1].fnuptake == Approx(1));

	// less N available than needed
	ncompete(indivs, 1.0);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));
	REQUIRE(indivs[1].fnuptake == Approx(0.5));

	// more N available than needed
	ncompete(indivs, 3.0);
	
	REQUIRE(indivs[0].fnuptake == Approx(1));
	REQUIRE(indivs[1].fnuptake == Approx(1));

	// less N available than needed
	ncompete(indivs, 1.0);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));
	REQUIRE(indivs[1].fnuptake == Approx(0.5));

	// make the second indiv twice as strong
	indivs[1].strength *= 2;

	// less N available than needed
	ncompete(indivs, 1.0);

	REQUIRE(indivs[0].fnuptake == Approx(1.0/3.0));
	REQUIRE(indivs[1].fnuptake == Approx(2.0/3.0));

	// reduce the demand for the stronger individual
	indivs[1].ndemand = 0.1;

	// less N available than needed
	ncompete(indivs, 1.0);

	REQUIRE(indivs[0].fnuptake == Approx(0.9));
	REQUIRE(indivs[1].fnuptake == Approx(1));

}

TEST_CASE("ncompete/triple", "Testing three individuals") {
	std::vector<NCompetingIndividual> indivs(3);

    indivs[0].ndemand = 1;
    indivs[0].strength = 1;

    indivs[1].ndemand = 0.01;
    indivs[1].strength = 5;

    indivs[2].ndemand = 5;
    indivs[2].strength = 50;

	ncompete(indivs, 5); 

	REQUIRE(indivs[0].fnuptake == Approx(0.097843137)); // was 0.24
	REQUIRE(indivs[1].fnuptake == Approx(1));
	REQUIRE(indivs[2].fnuptake == Approx(0.97843137));

    indivs[0].ndemand = 0.01;
    indivs[0].strength = 5;

    indivs[1].ndemand = 5;
    indivs[1].strength = 50;

	indivs[2].ndemand = 1;
    indivs[2].strength = 1;

    ncompete(indivs, 5); 

	REQUIRE(indivs[0].fnuptake == Approx(1));
	REQUIRE(indivs[1].fnuptake == Approx(0.97843137));
	REQUIRE(indivs[2].fnuptake == Approx(0.097843137));

    indivs[0].ndemand = 5;
    indivs[0].strength = 50;

	indivs[1].ndemand = 1;
    indivs[1].strength = 1;

    indivs[2].ndemand = 0.01;
    indivs[2].strength = 5;

    ncompete(indivs, 5); 

	REQUIRE(indivs[0].fnuptake == Approx(0.97843137));
	REQUIRE(indivs[1].fnuptake == Approx(0.097843137));
	REQUIRE(indivs[2].fnuptake == Approx(1));
}

TEST_CASE("ncompete/four", "Testing four individuals") {
	std::vector<NCompetingIndividual> indivs(4);

	indivs[0].ndemand = 0.5;
    indivs[0].strength = 2;

	indivs[1].ndemand = 8;
    indivs[1].strength = 5;

    indivs[2].ndemand = 0.5;
    indivs[2].strength = 0.1;

	indivs[3].ndemand = 2;
    indivs[3].strength = 100;

    ncompete(indivs, 10); 

	REQUIRE(indivs[0].fnuptake == Approx(1));
	REQUIRE(indivs[1].fnuptake == Approx(0.919118));
	REQUIRE(indivs[2].fnuptake == Approx(0.294118));
	REQUIRE(indivs[3].fnuptake == Approx(1));

    indivs[0].ndemand = 0.5;
    indivs[0].strength = 0.1;

	indivs[1].ndemand = 2;
    indivs[1].strength = 100;

	indivs[2].ndemand = 0.5;
    indivs[2].strength = 2;

	indivs[3].ndemand = 8;
    indivs[3].strength = 5;

    ncompete(indivs, 10); 

	REQUIRE(indivs[0].fnuptake == Approx(0.294118));
	REQUIRE(indivs[1].fnuptake == Approx(1));
	REQUIRE(indivs[2].fnuptake == Approx(1));
	REQUIRE(indivs[3].fnuptake == Approx(0.919118));
}
