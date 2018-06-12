#define BOOST_TEST_MODULE "DOCI_orbital_optimization_test"


#include <hf.hpp>
#include <cpputil.hpp>

#include "FCI.hpp"


#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



// dim = 2 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_sto_3g ) {

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
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());

    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.orbitalOptimize(&davidson_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;

    std::cout << "reference FCI energy: " << std::setprecision(15) << reference_fci_energy << std::endl;
    std::cout << "OO-DOCI energy: " << std::setprecision(15) << OO_DOCI_energy << std::endl;

    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}


// dim = 4 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31g ) {

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
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());

    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.orbitalOptimize(&davidson_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;

    std::cout << "reference FCI energy: " << std::setprecision(15) << reference_fci_energy << std::endl;
    std::cout << "OO-DOCI energy: " << std::setprecision(15) << OO_DOCI_energy << std::endl;


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
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());

    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.orbitalOptimize(&davidson_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;

    std::cout << "reference FCI energy: " << std::setprecision(15) << reference_fci_energy << std::endl;
    std::cout << "OO-DOCI energy: " << std::setprecision(15) << OO_DOCI_energy << std::endl;


    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( tessteat ) {


    libwint::Molecule beh2 ("../tests/beh2.xyz");
    libwint::AOBasis ao_basis (beh2, "6-31g**");
    ao_basis.calculateIntegrals();

    std::cout << "K: " << ao_basis.calculateNumberOfBasisFunctions() << std::endl;
}

// H4 (4e-)
// 6-31G**
//      K:      20
//      DOCI:   190
//      FCI:    36.100

// aug-ccpVDZ
//      K:      36
//      DOCI:   630
//      FCI:    396.900


// BeH2 (4e-)
// 6-31G**
//      K:      25
//      DOCI:   300
//      FCI:    90.000

// aug-ccpVDZ
//      K:      41
//      DOCI:   820
//      FCI:    672.400




// dim = 21
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2o_STO_3G ) {

    double reference_fci_energy = -82.9840637947819;

    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2o ("../tests/reference_data/h2o.xyz");
    double internuclear_repulsion_energy = h2o.calculateInternuclearRepulsionEnergy();  // 8.00236697416617

    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2o, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


//    // Get the FCI natural orbitals
//    ci::FCI fci (so_basis, 5, 5);  // N_alpha = 5, N_beta = 5
//    // Specify solver options and solve the eigenvalue problem
//    numopt::eigenproblem::DenseSolverOptions dense_options;
//    fci.solve(&dense_options);
//    std::cout << "FCI ENERGY DENSE CALCULATED : " << std::setprecision(15) << fci.get_lowest_eigenvalue() << std::endl;
//
//    fci.calculate1RDMs();
//    Eigen::MatrixXd D = fci.get_one_rdm();
//    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (D);
//
//    Eigen::MatrixXd U = saes.eigenvectors();
//    so_basis.rotate(U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    ci::DOCI doci (so_basis, h2o);

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());

    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.orbitalOptimize(&davidson_options);
    double OO_DOCI_energy = doci.get_lowest_eigenvalue() + internuclear_repulsion_energy;

//    std::cout << "reference FCI energy: " << std::setprecision(15) << fci.get_lowest_eigenvalue() << std::endl;
    std::cout << "OO-DOCI energy: " << std::setprecision(15) << doci.get_lowest_eigenvalue() << std::endl;


    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}
