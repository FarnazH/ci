#define BOOST_TEST_MODULE "DOCI_orbital_optimization_test"


#include <hf.hpp>
#include <cpputil.hpp>

#include "FCI.hpp"


#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



// dim = 2 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_sto_3g ) {

    // Check if OO-DOCI = FCI for a two-electron system
    double reference_fci_energy = -1.13726333769813;


    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327

    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Do the DOCI orbital optimization
    ci::DOCI doci (so_basis, h2);

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci.orbitalOptimize(&dense_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}


// dim = 4 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31g ) {

    // Check if OO-DOCI = FCI for a two-electron system, starting from the FCI naturals
    double reference_fci_energy = -1.15168629203274;


    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327

    libwint::AOBasis ao_basis (h2, "6-31G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Get the FCI natural orbitals
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions fci_dense_options;
    fci.solve(&fci_dense_options);

    fci.calculate1RDMs();
    Eigen::MatrixXd D = fci.get_one_rdm();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (D);

    Eigen::MatrixXd U = saes.eigenvectors();
    so_basis.rotate(U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    ci::DOCI doci (so_basis, h2);

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DenseSolverOptions doci_dense_options;
    doci.orbitalOptimize(&doci_dense_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}


// dim = 10 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31gxx ) {

    double reference_fci_energy = -1.16514875501195;

    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327

    libwint::AOBasis ao_basis (h2, "6-31G**");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Get the FCI natural orbitals
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

    fci.calculate1RDMs();
    Eigen::MatrixXd D = fci.get_one_rdm();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (D);

    Eigen::MatrixXd U = saes.eigenvectors();
    so_basis.rotate(U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    ci::DOCI doci (so_basis, h2);

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DenseSolverOptions doci_dense_options;
    doci.orbitalOptimize(&doci_dense_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}


// dim = 10 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31gxx_Davidson ) {

    double reference_fci_energy = -1.16514875501195;

    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327

    libwint::AOBasis ao_basis (h2, "6-31G**");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Get the FCI natural orbitals
    ci::FCI fci (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

    fci.calculate1RDMs();
    Eigen::MatrixXd D = fci.get_one_rdm();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (D);

    Eigen::MatrixXd U = saes.eigenvectors();
    so_basis.rotate(U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    ci::DOCI doci (so_basis, h2);

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());
    initial_guess(0) = 1;  // RHF initial guess
    davidson_options.X_0 = initial_guess;
    doci.orbitalOptimize(&davidson_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}
