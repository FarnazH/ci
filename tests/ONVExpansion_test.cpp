#define BOOST_TEST_MODULE "ONVExpansion_test"


#include "ONVExpansion.hpp"

#include <hf.hpp>


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( constructor_initializer_list ) {

    // Create a faulty expansion: one of the orbitals is different
    BOOST_CHECK_THROW((ci::ONVExpansion<unsigned long> {{bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (1, 3), 1.0},
                                                        {bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (2, 2), 0.0}}),
                      std::invalid_argument);


    // Create a correct expansion
    BOOST_CHECK_NO_THROW((ci::ONVExpansion<unsigned long> {{bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), 1.0},
                                                           {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.0}}));
}


BOOST_AUTO_TEST_CASE ( isEqual ) {

    ci::ONVExpansion<unsigned long> expansion1 {{bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (1, 3), 1.0},
                                                {bmqc::SpinString<unsigned long> (2, 3), bmqc::SpinString<unsigned long> (2, 3), 0.0}};

    ci::ONVExpansion<unsigned long> expansion2 {{bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (1, 3), 1.0},
                                                {bmqc::SpinString<unsigned long> (2, 3), bmqc::SpinString<unsigned long> (2, 3), 0.0},
                                                {bmqc::SpinString<unsigned long> (3, 3), bmqc::SpinString<unsigned long> (3, 3), 0.0}};

    ci::ONVExpansion<unsigned long> expansion3 {{bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), 1.0},
                                                {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.0}};

    ci::ONVExpansion<unsigned long> expansion4 {{bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (1, 3), 0.707},
                                                {bmqc::SpinString<unsigned long> (2, 3), bmqc::SpinString<unsigned long> (2, 3), 0.707}};

    ci::ONVExpansion<unsigned long> expansion5 {{bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (1, 3), -1.0},
                                                {bmqc::SpinString<unsigned long> (2, 3), bmqc::SpinString<unsigned long> (2, 3), 0.0}};

    BOOST_CHECK(!expansion1.isEqual(expansion2));  // wrong dimension
    BOOST_CHECK(!expansion1.isEqual(expansion3));  // different SpinStrings
    BOOST_CHECK(!expansion1.isEqual(expansion4));  // different coefficients

    BOOST_CHECK(expansion1.isEqual(expansion5));  // same eigenvector
}


BOOST_AUTO_TEST_CASE ( constructor_filename ) {

    // Create a reference expansion and a test expansion that is read in from a file
    ci::ONVExpansion<unsigned long> ref_expansion {{bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (1, 46), 1.0},
                                                   {bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (2, 46), 0.0}};

    ci::ONVExpansion<unsigned long> test_expansion ("../tests/reference_data/test_GAMESS_expansion");


    // Check if both expansions are considered equal
    BOOST_CHECK(ref_expansion.isEqual(test_expansion, 1.0e-06));
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
//
//
//BOOST_AUTO_TEST_CASE ( rdm1 ) {
//
//    // Do an H2@FCI//STO-3G calculation
//
//    // Prepare the AO basis
//    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
//    libwint::AOBasis ao_basis (h2, "STO-3G");
//    ao_basis.calculateIntegrals();
//
//    // Prepare the SO basis from RHF coefficients
//    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
//    rhf.solve();
//    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());
//
//
//    // Do a dense FCI calculation based on a given SO basis
//    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
//    fci.solve(numopt::eigenproblem::SolverType::DENSE);
//
//
//    fci.calculate1RDMs();
//    std::cout << fci.get_one_rdm_aa() << std::endl;
//
//
//    // Create a reference expansion
//    ci::ONVExpansion<unsigned long> ref_expansion {{bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), -0.993601},
//                                                   {bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (2, 2), -2.40468e-16},
//                                                   {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (1, 2), -3.04909e-16},
//                                                   {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.112949}};
//
//
//    ref_expansion.calculate1RDMs();
//}