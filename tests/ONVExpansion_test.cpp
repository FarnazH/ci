#define BOOST_TEST_MODULE "ONVExpansion_test"


#include "ONVExpansion.hpp"



#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#include <libint2/basis.h>


BOOST_AUTO_TEST_CASE ( filename_constructor ) {

    // Create a reference expansion and a test expansion that is read in from a file
    ci::ONVExpansion<unsigned long> ref_expansion = {{bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (1, 46), 1.0},
                                                     {bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (2, 46), 0.0}};

    ci::ONVExpansion<unsigned long> test_expansion ("../tests/reference_data/test_GAMESS_expansion");


    // Check if both expansions are considered equal
    for (size_t i = 0; i < ref_expansion.size(); i++) {
        BOOST_CHECK(test_expansion[i].isEqual(ref_expansion[i]));
    }
}
