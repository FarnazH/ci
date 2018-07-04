#define BOOST_TEST_MODULE "ONVExpansionComponent_test"


#include "ONVExpansionComponent.hpp"



#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( isEqual ) {

    bmqc::SpinString<unsigned long> a1 (1, 2);  // "01"
    bmqc::SpinString<unsigned long> b1 (1, 2);  // "01"
    ci::ONVExpansionComponent<unsigned long> component1 {a1, b1, 1.000};

    bmqc::SpinString<unsigned long> a2 (1, 2);  // "01"
    bmqc::SpinString<unsigned long> b2 (1, 2);  // "01"
    ci::ONVExpansionComponent<unsigned long> component2 {a2, b2, 1.000};

    bmqc::SpinString<unsigned long> a3 (1, 2);  // "01"
    bmqc::SpinString<unsigned long> b3 (1, 3);  // "001"
    ci::ONVExpansionComponent<unsigned long> component3 {a3, b3, 1.000};

    bmqc::SpinString<unsigned long> a4 (1, 2);  // "01"
    bmqc::SpinString<unsigned long> b4 (1, 2);  // "01"
    ci::ONVExpansionComponent<unsigned long> component4 {a4, b4, 0.500};


    BOOST_CHECK(component1.isEqual(component2));
    BOOST_CHECK(!component1.isEqual(component3));
    BOOST_CHECK(!component1.isEqual(component4));
}
