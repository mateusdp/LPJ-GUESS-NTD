///////////////////////////////////////////////////////////////////////////////////////
/// \file date_test.cpp
/// \brief Unit tests for CF time functionality
///
/// \author Joe Siltberg
/// $Date: 2019-11-06 17:45:44 +0100 (Wed, 06 Nov 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "cftime.h"
#include <time.h>

#ifdef HAVE_NETCDF

using namespace GuessNC::CF;

namespace {

bool verify_datetime(const DateTime& dt, 
                     int year, int month, int day, int hour, int minute, double second) {
	return dt.get_year() == year &&
		dt.get_month() == month &&
		dt.get_day() == day &&
		dt.get_hour() == hour &&
		dt.get_minute() == minute &&
		dt.get_second() == Approx(second);
}

}

TEST_CASE("CF::DateTime/construction", "Tests construction of DateTime objects") {
	// default constructor
	REQUIRE(verify_datetime(DateTime(), 
	                        0, 1, 1, 0, 0, 0));

	// string constructor
	REQUIRE(verify_datetime(DateTime("2013-01-25 09:09:09"),
	                        2013,01,25,9,9,9));

	// fractional seconds
	REQUIRE(verify_datetime(DateTime("1970-01-01 00:00:00.5"),
	                        1970,01,01,0,0,.5));

	// date only
	REQUIRE(verify_datetime(DateTime("1900-10-01"),
	                        1900,10,01,0,0,0));
}

TEST_CASE("CF::DateTime/add_time", "Tests the DateTime::add_time function") {

	// test NO_LEAP calendar

	DateTime dt("1980-01-01");
	
	dt.add_time(1, DAYS, NO_LEAP);
	REQUIRE(verify_datetime(dt, 1980,01,02,0,0,0));
	
	dt.add_time(30, DAYS, NO_LEAP);
	REQUIRE(verify_datetime(dt, 1980,02,01,0,0,0));
	
	dt.add_time(28, DAYS, NO_LEAP);
	REQUIRE(verify_datetime(dt, 1980,03,01,0,0,0));

	dt.add_time(365, DAYS, NO_LEAP);
	REQUIRE(verify_datetime(dt, 1981,03,01,0,0,0));

	dt.add_time(365*24+0.5, HOURS, NO_LEAP);
	REQUIRE(verify_datetime(dt, 1982,03,01,0,30,0));

	dt.add_time(0.1, SECONDS, NO_LEAP); // test rounding
	REQUIRE(verify_datetime(dt, 1982,03,01,0,30,0));

	dt.add_time(23.5*3600+365*24*3600*200.0, SECONDS, NO_LEAP);
	REQUIRE(verify_datetime(dt, 2182,03,02,0,0,0));

	// test STANDARD calendar
	dt = DateTime("1980-01-01");
	
	dt.add_time(1, DAYS, STANDARD);
	REQUIRE(verify_datetime(dt, 1980,01,02,0,0,0));

	dt.add_time(30, DAYS, STANDARD);
	REQUIRE(verify_datetime(dt, 1980,02,01,0,0,0));

	dt.add_time(28, DAYS, STANDARD);
	REQUIRE(verify_datetime(dt, 1980,02,29,0,0,0));

	dt.add_time(1, DAYS, STANDARD);
	REQUIRE(verify_datetime(dt, 1980,03,01,0,0,0));

	dt.add_time(364, DAYS, STANDARD);
	REQUIRE(verify_datetime(dt, 1981,02,28,0,0,0));

	dt.add_time(1, DAYS, STANDARD);
	REQUIRE(verify_datetime(dt, 1981,03,01,0,0,0));

	dt.add_time(366, DAYS, STANDARD);
	REQUIRE(verify_datetime(dt, 1982,03,02,0,0,0));

	// step 365*25 days forward and use C library time functions
	// (which use a gregorian calendar) to verify we reach the
	// correct date
	dt.add_time(365*25, DAYS, STANDARD);

	tm start_tm;
	start_tm.tm_sec = 0;
	start_tm.tm_min = 0;
	start_tm.tm_hour = 0;
	start_tm.tm_mday = 2;
	start_tm.tm_mon = 2;
	start_tm.tm_year = 82;
	start_tm.tm_isdst = 0;

	time_t end_t = mktime(&start_tm)+25*365*24*3600;
	tm* end_tm = localtime(&end_t);

	REQUIRE(verify_datetime(dt, 
	                        1900+end_tm->tm_year, end_tm->tm_mon+1, end_tm->tm_mday,
	                        end_tm->tm_hour, end_tm->tm_min, end_tm->tm_sec));
}

TEST_CASE("CF::TimeUnitSpecification", "Tests CF time unit specifications") {
	
	TimeUnitSpecification tus("days since 1970-01-01 00:00:00");
	
	REQUIRE(verify_datetime(tus.get_date_time(1, STANDARD),
	                        1970,01,02,0,0,0));

	tus = TimeUnitSpecification("hours since 1980-02-28 23:00:00");
	REQUIRE(verify_datetime(tus.get_date_time(1, STANDARD),
	                        1980,02,29,0,0,0));
}

#endif // HAVE_NETCDF
