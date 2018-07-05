#define BOOST_TEST_MODULE "ONVExpansion_test"


#include "ONVExpansion.hpp"

#include <hf.hpp>


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( constructor_filename ) {

    // Create a reference expansion and a test expansion that is read in from a file
    ci::ONVExpansion<unsigned long> ref_expansion = {{bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (1, 46), 1.0},
                                                     {bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (2, 46), 0.0}};

    ci::ONVExpansion<unsigned long> test_expansion ("../tests/reference_data/test_GAMESS_expansion");


    // Check if both expansions are considered equal
    for (size_t i = 0; i < ref_expansion.size(); i++) {
        BOOST_CHECK(test_expansion[i].isEqual(ref_expansion[i]));
    }
}


BOOST_AUTO_TEST_CASE ( constructor_fci ) {

    // Do an H2@FCI//STO-3G calculation

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    // Create the reference expansion (constructing the ONVExpansion manually)
    ci::ONVExpansion<unsigned long> ref_expansion {{bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), -0.993601},
                                                   {bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (2, 2), -2.40468e-16},
                                                   {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (1, 2), -3.04909e-16},
                                                   {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.112949}};


    // Create the test expansion
    ci::ONVExpansion<unsigned long> test_expansion (fci);


    BOOST_CHECK(test_expansion.isEqual(ref_expansion, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( constructor_doci ) {

    // Do an H2@DOCI//STO-3G calculation

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::DOCI doci (so_basis, 2);  // 2 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);


    // Set up the reference ONV expansion
    ci::ONVExpansion<unsigned long> ref_expansion {{bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), -0.993601},
                                                   {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.112949}};


    // Create the test expansion
    ci::ONVExpansion<unsigned long> test_expansion (doci);


    BOOST_CHECK(test_expansion.isEqual(ref_expansion, 1.0e-06));
}
