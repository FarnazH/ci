#define BOOST_TEST_MODULE "DOCI_test"

#include "DOCI.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( DOCI_co_klaas ) {

    // Klaas' reference DOCI energy for CO
    

}



BOOST_AUTO_TEST_CASE ( DOCI_h2o_klaas ) {

    // Klaas' reference DOCI energy for water


    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);
    doci::CI_basis ciBasis (rhf);

    doci::DOCI doci_test = doci::DOCI(&ciBasis);

    const doci::State& ground = doci_test.getLowestEigenState();

    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    BOOST_CHECK(std::abs(en - (-74.9771)) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE ( DOCI_beh_klaas ) {

    //beh_cation
    const std::string doci = "../tests/reference_data/doci/beh_cation_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::DOCI doci_test(&ciBasis);
    //getting the electronic ground energy
    const doci::State &ground = doci_test.getLowestEigenState();
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-14.8782216937)) < 1.0e-06); //energy from beh cation

}

BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas ) {
    //lih
    const std::string doci = "../tests/reference_data/doci/lih_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::DOCI doci_test(&ciBasis);
    //getting the electronic ground energy
    const doci::State &ground = doci_test.getLowestEigenState();
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-8.0029560313)) < 1.0e-06); //energy from lih
}
