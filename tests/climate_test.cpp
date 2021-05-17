///////////////////////////////////////////////////////////////////////////////////////
/// \file climate_test.cpp
/// \brief Unit tests for functions processing climate data
///
/// \author Joe Siltberg
/// $Date: 2014-03-03 08:50:44 +0100 (Mon, 03 Mar 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "catch.hpp"

#include "driver.h"
#include <algorithm>
#include <vector>

namespace {

/// Convenience function for testing prdaily
/** Tests prdaily for given monthly conditions, the same conditions
  * are used throughout the year, makes it easy to write test cases
  * without specifying 24 values.
  */
bool verify_prdaily_single_month(double prec, double wetdays) {

	double monthly_prec[12];
	std::fill(monthly_prec, monthly_prec+12, prec);

	double monthly_wetdays[12];
	std::fill(monthly_wetdays, monthly_wetdays+12, wetdays);

	double days[365];

	long seed = 12345678;
	prdaily(monthly_prec, days, monthly_wetdays, seed, 
			  false /* truncate set to false since we want to verify the sum */
			  );

	// Verify monthly sums and number of wet days
	const double SUM_TOLERANCE = 0.1;

	Date date;
	int current_day = 0;
	for (int m = 0; m < 12; m++) {
		double sum = 0;
		int wetcount = 0;

		for (int d = 0; d < date.ndaymonth[m]; d++) {
			sum += days[current_day];
			if (days[current_day] > 0) {
				wetcount++;
			}
			current_day++;
		}

		// Check that the sum for this month isn't too far from
		// the prescribed monthly precipitation
		if (fabs(sum-prec) > SUM_TOLERANCE) {
			return false;
		}

		// Verify wetcount?
	}

	return true;
}

bool verify_interp_monthly_means_conserve(const double* mvals,
                                          double minimum = -std::numeric_limits<double>::max(),
                                          double maximum = std::numeric_limits<double>::max()) {

	const double TOLERANCE = 0.0001;

	// Set upper and lower limits for allowed daily values
	// (apart from the limits supplied as parameters)
	std::vector<double> upper_limit(12), lower_limit(12);

	for (int m = 0; m < 12; m++) {
		int next_month = (m+1)%12;
		int prev_month = (m+11)%12;

		double next = mvals[next_month];
		double prev = mvals[prev_month];
		double current = mvals[m];

		double smallest = std::min(current, std::min(next, prev));
		double largest = std::max(current, std::max(next, prev));

		upper_limit[m] = std::max(largest, current+(current-smallest));
		lower_limit[m] = std::min(smallest, current-(largest-current));
	}

	double dvals[365];

	interp_monthly_means_conserve(mvals, dvals, minimum, maximum);

	std::vector<double> sums(12, 0);

	// Make sure daily values are within allowed limits

	Date date;
	date.init(1);

	for (int i = 0; i < 365; i++) {
		sums[date.month] += dvals[i];

		if (dvals[i] > upper_limit[date.month] + TOLERANCE ||
		    dvals[i] < lower_limit[date.month] - TOLERANCE) {
			return false;
		}

		if (dvals[i] < minimum || dvals[i] > maximum) {
			return false;
		}

		date.next();
	}

	// Make sure monthly means are conserved

	for (int m = 0; m < 12; m++) {
		if (fabs(sums[m]/date.ndaymonth[m] - mvals[m]) > TOLERANCE) {
			return false;
		}
	}

	return true;
}

}

TEST_CASE("climate/prdaily", "Tests for the prdaily function") {
	// Test no water
	REQUIRE(verify_prdaily_single_month(0, 0));

	// Regression test: 
	// Very small number of wet days used to cause infinite loop
	REQUIRE(verify_prdaily_single_month(1, 0.001));

	// Regression test: 
	// Very little precipitation and many wet days used to cause infinite loop
	REQUIRE(verify_prdaily_single_month(0.1, 30));
}

TEST_CASE("climate/interp_monthly_means_conserve", "Tests the monthly to daily interpolation") {

	double test1[] = { 0, 10, 20, 15, 15, 15, 40, 0, 40, 30, 20, 5};

	REQUIRE(verify_interp_monthly_means_conserve(test1, 0, 40));
}
