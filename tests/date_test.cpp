///////////////////////////////////////////////////////////////////////////////////////
/// \file date_test.cpp
/// \brief Unit tests for the Date class
///
/// \author Joe Siltberg
/// $Date: 2014-06-16 15:46:53 +0200 (Mon, 16 Jun 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "catch.hpp"

#include "guess.h"

namespace {

/// Help function that simply calls Date::next n times
void take_n_steps(Date& d, int n) {
	 for (int i = 0; i < n; i++) {
		  d.next();
	 }
}

}


TEST_CASE("date/construction", "Some basic tests of constructing dates") {
	 Date d;
	 d.init(1);

	 REQUIRE(d.day == 0);
	 REQUIRE(d.dayofmonth == 0);
	 REQUIRE(d.month == 0);
	 REQUIRE(d.year == 0);
	 REQUIRE(d.islastyear);
	 REQUIRE(!d.islastmonth);
	 REQUIRE(!d.islastday);
	 REQUIRE(!d.ismidday);
}


TEST_CASE("date/stepping", "Tests Date::next()") {
	 Date d;
	 d.init(1);

	 // Given that we take these number of steps...
	 int steps[] = { 1, 30, 28, 306 };
	 // ...we expect to end up at these combinations of year, month and day
	 int expectations[][3] = { {0, 0, 1}, { 0, 1, 0 }, {0, 2, 0}, {1, 0, 0} };

	 for (int i = 0; i < sizeof(steps)/sizeof(steps[0]); i++) {

		  take_n_steps(d, steps[i]);

		  REQUIRE(d.year == expectations[i][0]);
		  REQUIRE(d.month == expectations[i][1]);
		  REQUIRE(d.dayofmonth == expectations[i][2]);
	 }
}

TEST_CASE("date/leap", "Tests isleap") {
	REQUIRE(!Date::is_leap(1900));
	REQUIRE(!Date::is_leap(1975));
	REQUIRE(Date::is_leap(1904));
	REQUIRE(Date::is_leap(2000));
}

TEST_CASE("date/months", "Tests prevmonth and nextmonth") {
	 Date d;
	 d.init(1);

	 // We start in January
	 REQUIRE(d.nextmonth() == 1);
	 REQUIRE(d.prevmonth() == 11);

	 // Go to February
	 take_n_steps(d, 31);

	 REQUIRE(d.nextmonth() == 2);
	 REQUIRE(d.prevmonth() == 0);

	 // Go to December
	 take_n_steps(d, 330);

	 REQUIRE(d.nextmonth() == 0);
	 REQUIRE(d.prevmonth() == 10);

	 // Go to January the following year
	 take_n_steps(d, 10);

	 REQUIRE(d.nextmonth() == 1);
	 REQUIRE(d.prevmonth() == 11);
}

TEST_CASE("date/calendar_year", "Tests the calendar year concept") {
	Date d;
	d.init(1);

	// If calendar year isn't set, the simulation year and calendar year will be the same
	REQUIRE(d.get_calendar_year() == 0);

	d.set_first_calendar_year(1900);
	REQUIRE(d.get_calendar_year() == 1900);

	// Go to next year
	take_n_steps(d, 400);
	REQUIRE(d.year == 1);
	REQUIRE(d.get_calendar_year() == 1901);
}
