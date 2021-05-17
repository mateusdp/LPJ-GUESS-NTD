///////////////////////////////////////////////////////////////////////////////////////
/// \file string_test.cpp
/// \brief Unit tests functionality in guessstring.h
///
/// \author Joe Siltberg
/// $Date: 2014-02-13 10:31:34 +0100 (Thu, 13 Feb 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "catch.hpp"

#include "guessstring.h"

TEST_CASE("trim", "Tests for trim()") {
	REQUIRE(trim("") == "");

	REQUIRE(trim(" ") == "");
	
	REQUIRE(trim("   trim \t\n  ") == "trim");

	REQUIRE(trim("a b") == "a b");
	
	REQUIRE(trim(" a b ") == "a b");
}

TEST_CASE("to_upper", "Tests for to_upper()") {
	REQUIRE(to_upper("abc") == "ABC");

	REQUIRE(to_upper("a b c") == "A B C");

	REQUIRE(to_upper("1") == "1");
}

TEST_CASE("to_lower", "Tests for to_lower()") {
	REQUIRE(to_lower("AbC") == "abc");

	REQUIRE(to_lower("A B C") == "a b c");

	REQUIRE(to_lower("1") == "1");
}
