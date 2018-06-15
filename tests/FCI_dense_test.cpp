#define BOOST_TEST_MODULE "FCI_dense_test"

#include "FCI.hpp"
#include <RHF.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( test_random_rotation_diagonal_dense_fci ) {

    // Check if a random rotation has no effect on the sum of the diagonal elements

    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2o ("../tests/reference_data/h2o.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2o, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);
    size_t K = so_basis.get_K();


    // Specify solver options and do a DOCI calculation
    ci::FCI fci1 (so_basis, 5, 5);
    numopt::eigenproblem::DenseSolverOptions dense_options1;
    fci1.solve(&dense_options1);

    Eigen::VectorXd diagonal1 = fci1.get_diagonal();


    // Get a random unitary matrix by diagonalizing a random symmetric matrix
    Eigen::MatrixXd A_random = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd A_symmetric = A_random + A_random.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
    Eigen::MatrixXd U_random = unitary_solver.eigenvectors();


    // Rotate the SOBasis using the random unitary matrix
    so_basis.rotate(U_random);


    // Specify solver options and do a DOCI calculation
    ci::FCI fci2 (so_basis, 5 ,5);
    numopt::eigenproblem::DenseSolverOptions dense_options2;
    fci2.solve(&dense_options2);

    Eigen::VectorXd diagonal2 = fci2.get_diagonal();


    BOOST_CHECK(std::abs(diagonal1.sum() - diagonal2.sum()) < 1.0e-10);
}


// dim = 100
BOOST_AUTO_TEST_CASE ( FCI_H2_Cristina_dense ) {

    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_fci_energy = -1.1651486697;


    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    libwint::AOBasis ao_basis (h2, "6-31g**");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);


    // Calculate the total FCI energy
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}


// dim = 441
BOOST_AUTO_TEST_CASE ( FCI_H2O_Psi4_GAMESS_dense ) {

    // Psi4 and GAMESS' FCI energy
    double reference_fci_energy = -75.0129803939602;


    // Prepare the AO basis
    libwint::Molecule water ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 5, 5);  // N_alpha = 5, N_beta = 5
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

    // Calculate the total energy
    double internuclear_repulsion_energy = water.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}


// dim = 2116
BOOST_AUTO_TEST_CASE ( FCI_He_Cristina_dense ) {

    // Cristina's He FCI energy
    double reference_fci_energy = -2.902533599;


    // Prepare the AO basis
    libwint::Molecule helium ("../tests/reference_data/he_cristina.xyz");
    libwint::AOBasis ao_basis (helium, "aug-cc-pVQZ");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (helium, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);


    // Calculate (get) the total FCI energy
    double test_fci_energy = fci.get_eigenvalue();  // no internuclear repulsion energy for atoms

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}
