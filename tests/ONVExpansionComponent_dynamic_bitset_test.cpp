#define BOOST_TEST_MODULE "ONVExpansionComponent_dynamic_bitset_test"


#include "ONVExpansionComponent.hpp"



#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( isEqual_dynamic_bitset ) {

    bmqc::SpinString<boost::dynamic_bitset<>> a1 (boost::dynamic_bitset<> (std::string("01")));
    bmqc::SpinString<boost::dynamic_bitset<>> b1 (boost::dynamic_bitset<> (std::string("01")));
    ci::ONVExpansionComponent<boost::dynamic_bitset<>> component1 {a1, b1, 1.000};

    bmqc::SpinString<boost::dynamic_bitset<>> a2 (boost::dynamic_bitset<> (std::string("01")));
    bmqc::SpinString<boost::dynamic_bitset<>> b2 (boost::dynamic_bitset<> (std::string("01")));
    ci::ONVExpansionComponent<boost::dynamic_bitset<>> component2 {a2, b2, 1.000};

    bmqc::SpinString<boost::dynamic_bitset<>> a3 (boost::dynamic_bitset<> (std::string("01")));
    bmqc::SpinString<boost::dynamic_bitset<>> b3 (boost::dynamic_bitset<> (std::string("001")));
    ci::ONVExpansionComponent<boost::dynamic_bitset<>> component3 {a3, b3, 1.000};

    bmqc::SpinString<boost::dynamic_bitset<>> a4 (boost::dynamic_bitset<> (std::string("01")));
    bmqc::SpinString<boost::dynamic_bitset<>> b4 (boost::dynamic_bitset<> (std::string("01")));
    ci::ONVExpansionComponent<boost::dynamic_bitset<>> component4 {a4, b4, 0.500};


    BOOST_CHECK(component1.isEqual(component2));
    BOOST_CHECK(!component1.isEqual(component3));
    BOOST_CHECK(!component1.isEqual(component4));
}
