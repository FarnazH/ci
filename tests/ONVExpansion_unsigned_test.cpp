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


/*
 *  AUXILIARY FUNCTIONS
 */
/**
 *  @return the trace of the given 2-RDM @param: d(p,p,q,q)
 *
 *  In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
 *  Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves

 */
double trace_2RDM(const Eigen::Tensor<double, 4>& d) {
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

    auto K = static_cast<size_t>(d.dimension(1));

    double trace = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace += d(p,p,q,q);
        }
    }

    return trace;
}


double full_trace_2RDM(const Eigen::Tensor<double, 4>& d) {
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

    auto K = static_cast<size_t>(d.dimension(1));

    double trace = 0.0;
    for (size_t p = 0; p < K; p++) {
        trace += d(p,p,p,p);
    }

    return trace;
}


Eigen::MatrixXd reduce_2RDM(const Eigen::Tensor<double, 4>& d) {
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

    auto K = static_cast<size_t>(d.dimension(1));

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(K, K);
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                D(p,q) += d(p,q,r,r);
            }
        }
    }

    return D;
}


BOOST_AUTO_TEST_CASE ( cristina_rdms ) {

    // Initialize the matrices/tensors that are read in from test files
    Eigen::MatrixXd ref_one_rdm_aa = Eigen::MatrixXd::Zero(46,46);
    Eigen::MatrixXd ref_one_rdm_bb = Eigen::MatrixXd::Zero(46,46);
    Eigen::MatrixXd ref_one_rdm = Eigen::MatrixXd::Zero(46,46);


    Eigen::Tensor<double, 4> ref_two_rdm_aaaa (46,46,46,46);
    Eigen::Tensor<double, 4> ref_two_rdm_aabb (46,46,46,46);
    Eigen::Tensor<double, 4> ref_two_rdm_bbaa (46,46,46,46);
    Eigen::Tensor<double, 4> ref_two_rdm_bbbb (46,46,46,46);
    Eigen::Tensor<double, 4> ref_two_rdm (46,46,46,46);

    ref_two_rdm_aaaa.setZero();
    ref_two_rdm_aabb.setZero();
    ref_two_rdm_bbaa.setZero();
    ref_two_rdm_bbbb.setZero();
    ref_two_rdm.setZero();


    // Read in the 1- and 2-RDMs (time-consuming step)
    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm1s_aa", ref_one_rdm_aa);
    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm1s_bb", ref_one_rdm_bb);
    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm1", ref_one_rdm);

    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm2s_aaaa", ref_two_rdm_aaaa);
    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm2s_aabb", ref_two_rdm_aabb);
    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm2s_bbaa", ref_two_rdm_bbaa);
    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm2s_bbbb", ref_two_rdm_bbbb);
    cpputil::io::readArrayFromFile("../tests/reference_data/rdm_test/rdm2", ref_two_rdm);

    // The 2-RDMs outputted in the files have different array indices, so we will shuffle them accordingly
    // Note that .shuffle() returns a copy, and doesn't shuffle in-place
    Eigen::array<int, 4> shuffle_indices ({0, 2, 1, 3});  // ikjl -> ijkl
    Eigen::Tensor<double, 4> ref_two_rdm_aaaa_shuffled = ref_two_rdm_aaaa.shuffle(shuffle_indices);
    Eigen::Tensor<double, 4> ref_two_rdm_aabb_shuffled = ref_two_rdm_aabb.shuffle(shuffle_indices);
    Eigen::Tensor<double, 4> ref_two_rdm_bbaa_shuffled = ref_two_rdm_bbaa.shuffle(shuffle_indices);
    Eigen::Tensor<double, 4> ref_two_rdm_bbbb_shuffled = ref_two_rdm_bbbb.shuffle(shuffle_indices);
    Eigen::Tensor<double, 4> ref_two_rdm_shuffled = ref_two_rdm.shuffle(shuffle_indices);


    // Check the traces of the 1-RDMs that are read in
    BOOST_CHECK(std::abs(ref_one_rdm_aa.trace() - 1) < 1.0e-09);  // 1 alpha electron
    BOOST_CHECK(std::abs(ref_one_rdm_bb.trace() - 1) < 1.0e-09);  // 1 beta electron
    BOOST_CHECK(std::abs(ref_one_rdm.trace() - 2) < 1.0e-08);  // 2 electrons

    // Check the traces of the 2-RDMs that are read in
    BOOST_CHECK(std::abs(trace_2RDM(ref_two_rdm_aaaa_shuffled) - 0) < 1.0e-12);  // 0 pairs of alpha electrons (two_rdm_aaaa is all zeros anyways)
    BOOST_CHECK(std::abs(trace_2RDM(ref_two_rdm_aabb_shuffled) - 1) < 1.0e-09);  // 1 pair of alpha-beta electrons
    BOOST_CHECK(std::abs(trace_2RDM(ref_two_rdm_bbaa_shuffled) - 1) < 1.0e-09);  // 1 pair of alpha-beta electrons
    BOOST_CHECK(std::abs(trace_2RDM(ref_two_rdm_bbbb_shuffled) - 0) < 1.0e-12);  // 0 pairs of beta electrons (two_rdm_bbbb is all zeros anyways)
    BOOST_CHECK(std::abs(trace_2RDM(ref_two_rdm_shuffled) - 1) < 1.0e-09);  // 1 pair of electrons

    // Check if the read in 1-RDM can be calculated from the read in 2-RDM
    BOOST_CHECK(ref_one_rdm.isApprox(reduce_2RDM(2*ref_two_rdm_shuffled)));



    // Read in the reference expansion and calculate its 1- and 2-RDMs
    ci::ONVExpansion<unsigned long> expansion ("../tests/reference_data/rdm_test/he_asymp0.50_3_gamess");
    expansion.calculate1RDMs();
    expansion.calculate2RDMs();


    // Check the traces of the calculated 1-RDMs
    BOOST_CHECK(std::abs(expansion.get_one_rdm_aa().trace() - 1) < 1.0e-09);  // 1 alpha electron
    BOOST_CHECK(std::abs(expansion.get_one_rdm_bb().trace() - 1) < 1.0e-09);  // 1 beta electron
    BOOST_CHECK(std::abs(expansion.get_one_rdm().trace() - 2) < 1.0e-08);  // 2 electrons

    // Check the traces of the calculated 2-RDMs
    BOOST_CHECK(std::abs(trace_2RDM(expansion.get_two_rdm_aaaa()) - 0) < 1.0e-12);  // 0 pairs of alpha electrons (two_rdm_aaaa is all zeros anyways)
    BOOST_CHECK(std::abs(trace_2RDM(expansion.get_two_rdm_aabb()) - 1) < 1.0e-09);  // 1 pair of alpha-beta electrons
    BOOST_CHECK(std::abs(trace_2RDM(expansion.get_two_rdm_bbaa()) - 1) < 1.0e-09);  // 1 pair of alpha-beta electrons
    BOOST_CHECK(std::abs(trace_2RDM(expansion.get_two_rdm_bbbb()) - 0) < 1.0e-12);  // 0 pairs of beta electrons (two_rdm_bbbb is all zeros anyways)
    BOOST_CHECK(std::abs(trace_2RDM(expansion.get_two_rdm()) - 2) < 1.0e-08);  // 1 pair of electrons times 2

    // Check if the calculated 1-RDM can be calculated from the calculated 2-RDM
    BOOST_CHECK(expansion.get_one_rdm().isApprox(reduce_2RDM(expansion.get_two_rdm())));



    // Check if the calculated and read in 1-RDMs are equal
    BOOST_CHECK(ref_one_rdm_aa.isApprox(expansion.get_one_rdm_aa(), 1.0e-06));
    BOOST_CHECK(ref_one_rdm_bb.isApprox(expansion.get_one_rdm_bb(), 1.0e-06));
    BOOST_CHECK(ref_one_rdm.isApprox(expansion.get_one_rdm(), 1.0e-06));

    // Check if the calculated and read in 2-RDMs are equal
    BOOST_CHECK(cpputil::linalg::areEqual(ref_two_rdm_aaaa_shuffled, expansion.get_two_rdm_aaaa(), 1.0e-12));  // are all zero anyways
    BOOST_CHECK(cpputil::linalg::areEqual(ref_two_rdm_aabb_shuffled, expansion.get_two_rdm_aabb(), 1.0e-04));
    BOOST_CHECK(cpputil::linalg::areEqual(ref_two_rdm_bbaa_shuffled, expansion.get_two_rdm_bbaa(), 1.0e-04));
    BOOST_CHECK(cpputil::linalg::areEqual(ref_two_rdm_bbbb_shuffled, expansion.get_two_rdm_bbbb(), 1.0e-12));  // are all zero anyways
    BOOST_CHECK(cpputil::linalg::areEqual(2*ref_two_rdm_shuffled, expansion.get_two_rdm(), 1.0e-04));  // factor half difference in other full 2-RDM


    // Calculate energies based on contractions of the integrals and density matrices
    libwint::SOBasis so_basis ("../tests/reference_data/rdm_test/he_asymp0.50_3.FCIDUMP", 46);
    // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
    //      0.5 g(p q r s) d(p q r s)
    Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};


    // Calculate the energy based off the read in density matrices
    double energy_by_contraction_read = (so_basis.get_h_SO() * ref_one_rdm).trace();
    Eigen::Tensor<double, 0> contraction_read = 0.5 * (so_basis.get_g_SO()).contract(2*ref_two_rdm_shuffled, contractions);
    energy_by_contraction_read += contraction_read(0);  // as the contraction is a scalar (a tensor of rank 0), we should access by (0)
    BOOST_CHECK(std::abs(energy_by_contraction_read - -3.0203317600) < 1.0e-08);

    // Calculate the energy based off the calculated density matrices
    double energy_by_contraction_calculated = (so_basis.get_h_SO() * expansion.get_one_rdm()).trace();
    Eigen::Tensor<double, 0> contraction_calculated = 0.5 * (so_basis.get_g_SO()).contract(expansion.get_two_rdm(), contractions);
    energy_by_contraction_calculated += contraction_calculated(0);  // as the contraction is a scalar (a tensor of rank 0), we should access by (0)
    BOOST_CHECK(std::abs(energy_by_contraction_calculated - -3.0203317600) < 1.0e-09);
}
