#define BOOST_TEST_MODULE "ONVExpansion_unsigned_test"


#include "ONVExpansion.hpp"

#include <hf.hpp>
#include <cpputil.hpp>


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( constructor_initializer_list_unsigned ) {

    // Create a faulty expansion: one of the orbitals is different
    BOOST_CHECK_THROW((ci::ONVExpansion<unsigned long> {{bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (1, 3), 1.0},
                                                        {bmqc::SpinString<unsigned long> (1, 3), bmqc::SpinString<unsigned long> (2, 2), 0.0}}),
                      std::invalid_argument);


    // Create a correct expansion
    BOOST_CHECK_NO_THROW((ci::ONVExpansion<unsigned long> {{bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), 1.0},
                                                           {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.0}}));
}


BOOST_AUTO_TEST_CASE ( isEqual_unsigned ) {

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


BOOST_AUTO_TEST_CASE ( constructor_filename_unsigned ) {

    // Create a reference expansion and a test expansion that is read in from a file
    ci::ONVExpansion<unsigned long> ref_expansion {{bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (1, 46), 1.0},
                                                   {bmqc::SpinString<unsigned long> (1, 46), bmqc::SpinString<unsigned long> (2, 46), 0.0}};

    ci::ONVExpansion<unsigned long> test_expansion ("../tests/reference_data/test_GAMESS_expansion");


    // Check if both expansions are considered equal
    BOOST_CHECK(ref_expansion.isEqual(test_expansion, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( constructor_fci_two_electrons_unsigned ) {

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
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

    // Create the reference expansion (constructing the ONVExpansion manually)
    ci::ONVExpansion<unsigned long> ref_expansion {
        {bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), -0.993601},
        {bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (2, 2), -2.40468e-16},
        {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (1, 2), -3.04909e-16},
        {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.112949}};


    // Create the test expansion
    ci::ONVExpansion<unsigned long> test_expansion (fci);


    BOOST_CHECK(test_expansion.isEqual(ref_expansion, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( constructor_fci_four_electrons_unsigned ) {

    // Do an H2(2-)@FCI//6-31G calculation

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    libwint::AOBasis ao_basis (h2, "6-31G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 2, 2);  // N_alpha = 2, N_beta = 2
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);


    // Create the reference expansion (constructing the ONVExpansion manually from the FCI calculation)
    ci::ONVExpansion<unsigned long> ref_expansion {
        {bmqc::SpinString<unsigned long> (3, 4), bmqc::SpinString<unsigned long> (3, 4), -0.984527},
        {bmqc::SpinString<unsigned long> (3, 4), bmqc::SpinString<unsigned long> (5, 4), -1.5665e-15},
        {bmqc::SpinString<unsigned long> (3, 4), bmqc::SpinString<unsigned long> (6, 4), 0.103568},
        {bmqc::SpinString<unsigned long> (3, 4), bmqc::SpinString<unsigned long> (9, 4), -0.0267219},
        {bmqc::SpinString<unsigned long> (3, 4), bmqc::SpinString<unsigned long> (10, 4), 2.65415e-17},
        {bmqc::SpinString<unsigned long> (3, 4), bmqc::SpinString<unsigned long> (12, 4), 0.0166314},

        {bmqc::SpinString<unsigned long> (5, 4), bmqc::SpinString<unsigned long> (3, 4), -1.19959e-15},
        {bmqc::SpinString<unsigned long> (5, 4), bmqc::SpinString<unsigned long> (5, 4), 0.0274075},
        {bmqc::SpinString<unsigned long> (5, 4), bmqc::SpinString<unsigned long> (6, 4), 5.68127e-17},
        {bmqc::SpinString<unsigned long> (5, 4), bmqc::SpinString<unsigned long> (9, 4), 5.33411e-17},
        {bmqc::SpinString<unsigned long> (5, 4), bmqc::SpinString<unsigned long> (10, 4), -0.0147791},
        {bmqc::SpinString<unsigned long> (5, 4), bmqc::SpinString<unsigned long> (12, 4), 2.61811e-17},

        {bmqc::SpinString<unsigned long> (6, 4), bmqc::SpinString<unsigned long> (3, 4), 0.103568},
        {bmqc::SpinString<unsigned long> (6, 4), bmqc::SpinString<unsigned long> (5, 4), 1.85374e-16},
        {bmqc::SpinString<unsigned long> (6, 4), bmqc::SpinString<unsigned long> (6, 4), 0.046555},
        {bmqc::SpinString<unsigned long> (6, 4), bmqc::SpinString<unsigned long> (9, 4), -0.0314105},
        {bmqc::SpinString<unsigned long> (6, 4), bmqc::SpinString<unsigned long> (10, 4), 6.51894e-17},
        {bmqc::SpinString<unsigned long> (6, 4), bmqc::SpinString<unsigned long> (12, 4), -0.00771233},

        {bmqc::SpinString<unsigned long> (9, 4), bmqc::SpinString<unsigned long> (3, 4), -0.0267219},
        {bmqc::SpinString<unsigned long> (9, 4), bmqc::SpinString<unsigned long> (5, 4), -1.64427e-16},
        {bmqc::SpinString<unsigned long> (9, 4), bmqc::SpinString<unsigned long> (6, 4), -0.0314105},
        {bmqc::SpinString<unsigned long> (9, 4), bmqc::SpinString<unsigned long> (9, 4), 0.0202616},
        {bmqc::SpinString<unsigned long> (9, 4), bmqc::SpinString<unsigned long> (10, 4), 7.20922e-18},
        {bmqc::SpinString<unsigned long> (9, 4), bmqc::SpinString<unsigned long> (12, 4), 0.00564933},

        {bmqc::SpinString<unsigned long> (10, 4), bmqc::SpinString<unsigned long> (3, 4), -1.60053e-16},
        {bmqc::SpinString<unsigned long> (10, 4), bmqc::SpinString<unsigned long> (5, 4), -0.0147791},
        {bmqc::SpinString<unsigned long> (10, 4), bmqc::SpinString<unsigned long> (6, 4), 1.11285e-16},
        {bmqc::SpinString<unsigned long> (10, 4), bmqc::SpinString<unsigned long> (9, 4), 8.69277e-18},
        {bmqc::SpinString<unsigned long> (10, 4), bmqc::SpinString<unsigned long> (10, 4), 0.036638},
        {bmqc::SpinString<unsigned long> (10, 4), bmqc::SpinString<unsigned long> (12, 4), 6.48637e-17},

        {bmqc::SpinString<unsigned long> (12, 4), bmqc::SpinString<unsigned long> (3, 4), 0.0166314},
        {bmqc::SpinString<unsigned long> (12, 4), bmqc::SpinString<unsigned long> (5, 4), -7.09468e-19},
        {bmqc::SpinString<unsigned long> (12, 4), bmqc::SpinString<unsigned long> (6, 4), -0.00771233},
        {bmqc::SpinString<unsigned long> (12, 4), bmqc::SpinString<unsigned long> (9, 4), 0.00564933},
        {bmqc::SpinString<unsigned long> (12, 4), bmqc::SpinString<unsigned long> (10, 4), 3.89442e-17},
        {bmqc::SpinString<unsigned long> (12, 4), bmqc::SpinString<unsigned long> (12, 4), -0.00291466}
    };


    // Create the test expansion
    ci::ONVExpansion<unsigned long> test_expansion (fci);


    BOOST_CHECK(test_expansion.isEqual(ref_expansion, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( constructor_doci_two_electrons_unsigned ) {

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
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci.solve(&dense_options);


    // Set up the reference ONV expansion
    ci::ONVExpansion<unsigned long> ref_expansion {{bmqc::SpinString<unsigned long> (1, 2), bmqc::SpinString<unsigned long> (1, 2), -0.993601},
                                                   {bmqc::SpinString<unsigned long> (2, 2), bmqc::SpinString<unsigned long> (2, 2), 0.112949}};


    // Create the test expansion
    ci::ONVExpansion<unsigned long> test_expansion (doci);


    BOOST_CHECK(test_expansion.isEqual(ref_expansion, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( constructor_doci_four_electrons_unsigned ) {

    // Do an H2(2-)@DOCI//6-31G calculation

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    libwint::AOBasis ao_basis (h2, "6-31G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci.solve(&dense_options);


    // Create the reference expansion (constructing the ONVExpansion manually from the DOCI calculation)
    ci::ONVExpansion<unsigned long> ref_expansion {
        {bmqc::SpinString<unsigned long> (3, 4), bmqc::SpinString<unsigned long> (3, 4), -0.997032},
        {bmqc::SpinString<unsigned long> (5, 4), bmqc::SpinString<unsigned long> (5, 4), 0.0246847},
        {bmqc::SpinString<unsigned long> (6, 4), bmqc::SpinString<unsigned long> (6, 4), 0.0553457},
        {bmqc::SpinString<unsigned long> (9, 4), bmqc::SpinString<unsigned long> (9, 4), 0.0227064},
        {bmqc::SpinString<unsigned long> (10, 4), bmqc::SpinString<unsigned long> (10, 4), 0.0416207},
        {bmqc::SpinString<unsigned long> (12, 4), bmqc::SpinString<unsigned long> (12, 4), -0.00263289}
    };


    // Create the test expansion
    ci::ONVExpansion<unsigned long> test_expansion (doci);


    BOOST_CHECK(test_expansion.isEqual(ref_expansion, 1.0e-06));
}


// K = 4, dim = 16
BOOST_AUTO_TEST_CASE ( one_rdms_fci_H2_6_31G_unsigned ) {

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
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

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
BOOST_AUTO_TEST_CASE ( two_rdms_fci_H2_6_31G_unsigned ) {

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
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

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
BOOST_AUTO_TEST_CASE ( one_rdms_fci_H2O_STO_3G_unsigned ) {

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
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);


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
BOOST_AUTO_TEST_CASE ( two_rdms_fci_H2O_STO_3G_unsigned ) {

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
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);


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
