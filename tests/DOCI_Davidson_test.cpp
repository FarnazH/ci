#define BOOST_TEST_MODULE "DOCI_Davidson_test"


#include <hf.hpp>
#include <cpputil.hpp>

#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( DOCI_li2_klaas_check_matvec_dense_Davidson ) {

    // Check if the matrix-vector product from Davidson and from the full dense Hamiltonian gives the same answer
    // I used a vector equal to ones
    Eigen::VectorXd matvec_dense = Eigen::VectorXd::Zero(816);
    Eigen::VectorXd matvec_davidson = Eigen::VectorXd::Zero(816);

    cpputil::io::readVectorFromFile("../tests/reference_data/li2_321g_klaas_matvec_with_ones_dense.data", matvec_dense);
    cpputil::io::readVectorFromFile("../tests/reference_data/li2_321g_klaas_matvec_with_ones_davidson.data", matvec_davidson);


    BOOST_CHECK(matvec_dense.isApprox(matvec_davidson, 1.0e-08));
}


// dim = 2
BOOST_AUTO_TEST_CASE ( DOCI_h2_sto3g_dense_vs_Davidson ) {

    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();

    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-08);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Calculate the DOCI energy using full diagonalization of the dense Hamiltonian
    ci::DOCI doci_dense (so_basis, h2);

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci_dense.solve(&dense_options);
    double doci_energy_dense = doci_dense.get_eigenvalue() + internuclear_repulsion_energy;


    // Calculate the DOCI energy using the Davidson algorithm
    ci::DOCI doci_davidson (so_basis, h2);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci_davidson.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci_davidson.solve(&davidson_options);
    double doci_energy_davidson = doci_davidson.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(doci_energy_dense - doci_energy_davidson) < 1.0e-08);
}


// dim = 4
BOOST_AUTO_TEST_CASE ( DOCI_h2_631g_dense_vs_Davidson ) {

    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();

    libwint::AOBasis ao_basis (h2, "6-31G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-08);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Calculate the DOCI energy using full diagonalization of the dense Hamiltonian
    ci::DOCI doci_dense (so_basis, h2);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci_dense.solve(&dense_options);
    double doci_energy_dense = doci_dense.get_eigenvalue() + internuclear_repulsion_energy;


    // Calculate the DOCI energy using the Davidson algorithm
    ci::DOCI doci_davidson (so_basis, h2);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci_davidson.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci_davidson.solve(&davidson_options);
    double doci_energy_davidson = doci_davidson.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(doci_energy_dense - doci_energy_davidson) < 1.0e-08);
}


// dim = 21
BOOST_AUTO_TEST_CASE ( DOCI_h2o_sto3g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for H2O@STO-3G
    double reference_doci_energy = -74.9671366903;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/h2o_sto3g_klaas.FCIDUMP", 7);  // 7 SOs
    ci::DOCI doci (so_basis, 10);  // 10 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 9.7794061444134091E+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 120
BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for BeH+
    double reference_doci_energy = -14.8782216937;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 120
BOOST_AUTO_TEST_CASE ( DOCI_n2_sto3g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for N2
    double reference_doci_energy = -107.5813316864;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/n2_sto-3g_klaas.FCIDUMP", 10);  // 10 SOs
    ci::DOCI doci (so_basis, 14);  // 10 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 2.3786407766990290E+01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 120
BOOST_AUTO_TEST_CASE ( DOCI_lih_631g_klaas_Davidson ) {
    
    // Klaas' reference DOCI energy for LiH
    double reference_doci_energy = -8.0029560313;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 816
BOOST_AUTO_TEST_CASE ( DOCI_li2_321g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for Li2
    double reference_doci_energy = -15.1153976060;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/li2_321g_klaas.FCIDUMP", 18);  // 18 SOs
    ci::DOCI doci (so_basis, 6);  // 6 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 3.0036546888874875e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 1287
BOOST_AUTO_TEST_CASE ( DOCI_h2o_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for H2O
    double reference_doci_energy = -76.0125161011;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/h2o_631g_klaas.FCIDUMP", 13);  // 13 SOs
    ci::DOCI doci (so_basis, 10);  // 10 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 9.7794061444134091E+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 376740
BOOST_AUTO_TEST_CASE ( DOCI_lif_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for LiF
    double reference_doci_energy = -107.0007150075;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/lif_631g_klaas.FCIDUMP", 28);  // 28 SOs
    ci::DOCI doci (so_basis, 12);  // 12 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 9.1249103487674024e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 1184040
BOOST_AUTO_TEST_CASE ( DOCI_co_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for CO
    double reference_doci_energy = -112.8392190587;


    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/co_631g_klaas.FCIDUMP", 28);  // 28 SOs
    ci::DOCI doci (so_basis, 14);  // 14 electrons

    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.solve(&davidson_options);


    // Calculate the total energy
    double internuclear_repulsion_energy = 2.2141305786610879e+01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-8);
}
