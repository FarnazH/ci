#define BOOST_TEST_MODULE "ONVExpansion_test"


#include "ONVExpansion.hpp"

#include <hf.hpp>
#include <cpputil.hpp>


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


// K = 4, dim = 16
BOOST_AUTO_TEST_CASE ( one_rdms_fci_H2_6_31G ) {

    // Do an H2@FCI//6-31G calculation

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    libwint::AOBasis ao_basis (h2, "6-31G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    fci.calculate1RDMs();
    Eigen::MatrixXd ref_one_rdm_aa = fci.get_one_rdm_aa();
    Eigen::MatrixXd ref_one_rdm_bb = fci.get_one_rdm_bb();
    Eigen::MatrixXd ref_one_rdm = fci.get_one_rdm();


    // Read in the FCI expansion into an ONVExpansion, and calculate the 1-RDMs
    ci::ONVExpansion<unsigned long> expansion (fci);
    expansion.calculate1RDMs();

    Eigen::MatrixXd test_one_rdm_aa = expansion.get_one_rdm_aa();
    Eigen::MatrixXd test_one_rdm_bb = expansion.get_one_rdm_bb();
    Eigen::MatrixXd test_one_rdm = expansion.get_one_rdm();


    BOOST_CHECK(test_one_rdm_aa.isApprox(ref_one_rdm_aa));
    BOOST_CHECK(test_one_rdm_bb.isApprox(ref_one_rdm_bb));
    BOOST_CHECK(test_one_rdm.isApprox(ref_one_rdm));
}


// K = 4, dim = 16
BOOST_AUTO_TEST_CASE ( two_rdms_fci_H2_6_31G ) {

    // Do an H2@FCI//STO-3G calculation

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    libwint::AOBasis ao_basis (h2, "6-31G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    fci.calculate2RDMs();
    Eigen::Tensor<double, 4> ref_two_rdm_aaaa = fci.get_two_rdm_aaaa();
    Eigen::Tensor<double, 4> ref_two_rdm_aabb = fci.get_two_rdm_aabb();
    Eigen::Tensor<double, 4> ref_two_rdm_bbaa = fci.get_two_rdm_bbaa();
    Eigen::Tensor<double, 4> ref_two_rdm_bbbb = fci.get_two_rdm_bbbb();
    Eigen::Tensor<double, 4> ref_two_rdm = fci.get_two_rdm();


    // Read in the FCI expansion into an ONVExpansion, and calculate the 2-RDMs
    ci::ONVExpansion<unsigned long> expansion (fci);

    expansion.calculate2RDMs();
    Eigen::Tensor<double, 4> test_two_rdm_aaaa = expansion.get_two_rdm_aaaa();
    Eigen::Tensor<double, 4> test_two_rdm_aabb = expansion.get_two_rdm_aabb();
    Eigen::Tensor<double, 4> test_two_rdm_bbaa = expansion.get_two_rdm_bbaa();
    Eigen::Tensor<double, 4> test_two_rdm_bbbb = expansion.get_two_rdm_bbbb();
    Eigen::Tensor<double, 4> test_two_rdm = expansion.get_two_rdm();


    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_aaaa, ref_two_rdm_aaaa, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_aabb, ref_two_rdm_aabb, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_bbaa, ref_two_rdm_bbaa, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_bbbb, ref_two_rdm_bbbb, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm, ref_two_rdm, 1.0e-06));
}


// K = 7, dim = 441
BOOST_AUTO_TEST_CASE ( one_rdms_fci_H2O_STO_3G ) {

    // Do an H2O@FCI//STO-3G calculation

    // Prepare the AO basis
    libwint::Molecule h2o ("../tests/reference_data/h2o.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2o, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 5, 5);  // N_alpha = 5, N_beta = 5
    fci.solve(numopt::eigenproblem::SolverType::DENSE);


    fci.calculate1RDMs();
    Eigen::MatrixXd ref_one_rdm_aa = fci.get_one_rdm_aa();
    Eigen::MatrixXd ref_one_rdm_bb = fci.get_one_rdm_bb();
    Eigen::MatrixXd ref_one_rdm = fci.get_one_rdm();


    // Read in the FCI expansion into an ONVExpansion, and calculate the 1-RDMs
    ci::ONVExpansion<unsigned long> expansion (fci);
    expansion.calculate1RDMs();

    Eigen::MatrixXd test_one_rdm_aa = expansion.get_one_rdm_aa();
    Eigen::MatrixXd test_one_rdm_bb = expansion.get_one_rdm_bb();
    Eigen::MatrixXd test_one_rdm = expansion.get_one_rdm();


    BOOST_CHECK(test_one_rdm_aa.isApprox(ref_one_rdm_aa));
    BOOST_CHECK(test_one_rdm_bb.isApprox(ref_one_rdm_bb));
    BOOST_CHECK(test_one_rdm.isApprox(ref_one_rdm));
}


// K = 7, dim = 441
BOOST_AUTO_TEST_CASE ( two_rdms_fci_H2O_STO_3G ) {

    // Do an H2O@FCI//STO-3G calculation

    // Prepare the AO basis
    libwint::Molecule h2o ("../tests/reference_data/h2o.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2o, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 5, 5);  // N_alpha = 5, N_beta = 5
    fci.solve(numopt::eigenproblem::SolverType::DENSE);


    fci.calculate2RDMs();
    Eigen::Tensor<double, 4> ref_two_rdm_aaaa = fci.get_two_rdm_aaaa();
    Eigen::Tensor<double, 4> ref_two_rdm_aabb = fci.get_two_rdm_aabb();
    Eigen::Tensor<double, 4> ref_two_rdm_bbaa = fci.get_two_rdm_bbaa();
    Eigen::Tensor<double, 4> ref_two_rdm_bbbb = fci.get_two_rdm_bbbb();
    Eigen::Tensor<double, 4> ref_two_rdm = fci.get_two_rdm();


    // Read in the FCI expansion into an ONVExpansion, and calculate the 2-RDMs
    ci::ONVExpansion<unsigned long> expansion (fci);

    expansion.calculate2RDMs();
    Eigen::Tensor<double, 4> test_two_rdm_aaaa = expansion.get_two_rdm_aaaa();
    Eigen::Tensor<double, 4> test_two_rdm_aabb = expansion.get_two_rdm_aabb();
    Eigen::Tensor<double, 4> test_two_rdm_bbaa = expansion.get_two_rdm_bbaa();
    Eigen::Tensor<double, 4> test_two_rdm_bbbb = expansion.get_two_rdm_bbbb();
    Eigen::Tensor<double, 4> test_two_rdm = expansion.get_two_rdm();


    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_aaaa, ref_two_rdm_aaaa, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_aabb, ref_two_rdm_aabb, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_bbaa, ref_two_rdm_bbaa, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm_bbbb, ref_two_rdm_bbbb, 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(test_two_rdm, ref_two_rdm, 1.0e-06));
}
