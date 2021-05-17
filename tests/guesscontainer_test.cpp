///////////////////////////////////////////////////////////////////////////////////////
/// \file guesscontainer_test.cpp
/// \brief Unit tests GuessContainer
///
/// \author Joe Siltberg
/// $Date: 2014-03-19 09:46:54 +0100 (Wed, 19 Mar 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "catch.hpp"

#include "guesscontainer.h"

TEST_CASE("indexing", "Tests random access") {
	GuessContainer<int> container;

	container.push_back(new int(1));
	container.push_back(new int(2));

	REQUIRE(container[0] == 1);
	REQUIRE(container[1] == 2);
	REQUIRE(container.size() == 2);
}

class InstanceCountingClass {
public:
	static int living_instances;

	InstanceCountingClass() {
		++living_instances;
	}

	~InstanceCountingClass() {
		--living_instances;
	}
};

int InstanceCountingClass::living_instances = 0;

TEST_CASE("memory", "Tests memory management in GuessContainer") {

	GuessContainer<InstanceCountingClass> container;

	container.push_back(new InstanceCountingClass());
	container.push_back(new InstanceCountingClass());
	container.push_back(new InstanceCountingClass());
	
	REQUIRE(InstanceCountingClass::living_instances == 3);

	container.erase(container.begin());

	REQUIRE(InstanceCountingClass::living_instances == 2);

	container.clear();

	REQUIRE(InstanceCountingClass::living_instances == 0);
}

TEST_CASE("iteration", "Tests iteration in GuessContainer") {
	GuessContainer<int> container;

	container.push_back(new int(1));
	container.push_back(new int(2));

	GuessContainer<int>::iterator itr = container.begin();

	REQUIRE(itr == container.begin());
	REQUIRE(itr++ == container.begin());
	REQUIRE(*itr == 2);

	itr = container.begin();

	REQUIRE(++itr != container.begin());
	REQUIRE(*itr == 2);
	
	REQUIRE(++itr == container.end());
}

TEST_CASE("erase", "Tests GuessContainer<T>::erase()") {
	GuessContainer<int> container;

	container.push_back(new int(1));
	container.push_back(new int(2));
	container.push_back(new int(3));
	
	GuessContainer<int>::iterator itr = container.erase(container.begin());

	REQUIRE(itr == container.begin());

	++itr;

	itr = container.erase(itr);

	REQUIRE(itr == container.end());
	
	itr = container.begin();
	itr = container.erase(itr);
	
	REQUIRE(itr == container.end());
	REQUIRE(container.size() == 0);
}
