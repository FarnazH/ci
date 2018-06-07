#define BOOST_TEST_MODULE "DOCI_orbital_optimization_test"


#include <hf.hpp>
#include <cpputil.hpp>

#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



// H2 FCI
// 6-31g**: -1.16514875501195
// 6-31g: -1.15168629203274
// STO-3G: -1.13726333769813


// dim = 2
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


    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-12);
}


// dim = 10
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_cristina ) {

    // Cristina's reference DOCI energy for H2
    // These are actually FCI results, but OO-DOCI for a 2-electron system should be exact within the given basis
    double reference_OO_DOCI_energy = -1.1651486697;  // -1,1697932197


    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327

    libwint::AOBasis ao_basis (h2, "6-31g**");

    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);

    ci::DOCI doci (so_basis, h2);

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());


    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;


    doci.orbitalOptimize(&davidson_options);






    // BOOST_CHECK(std::abs(test_doci_energy - reference_OO_DOCI_energy) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( OO_DOCI_h2o ) {

    // Cristina's reference DOCI energy for H2
    // These are actually FCI results, but OO-DOCI for a 2-electron system should be exact within the given basis
//    double reference_OO_DOCI_energy = -1.1651486697;  // -1,1697932197


    // Prepare an SOBasis from an RHF calculation
    libwint::Molecule h2o ("../tests/reference_data/h2o.xyz");
    double internuclear_repulsion_energy = h2o.calculateInternuclearRepulsionEnergy();

    libwint::AOBasis ao_basis (h2o, "STO-3G");

    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2o, ao_basis, 1.0e-06);
    rhf.solve();

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);

    ci::DOCI doci (so_basis, h2o);

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());

//
//    initial_guess(0) = 1;
//    davidson_options.X_0 = initial_guess;


    doci.orbitalOptimize(&davidson_options);




    // BOOST_CHECK(std::abs(test_doci_energy - reference_OO_DOCI_energy) < 1.0e-06);
}



BOOST_AUTO_TEST_CASE ( OO_DOCI_beh_cation ) {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons

    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());


    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;


    doci.orbitalOptimize(&davidson_options);



    // BOOST_CHECK(std::abs(test_doci_energy - reference_OO_DOCI_energy) < 1.0e-06);
}
