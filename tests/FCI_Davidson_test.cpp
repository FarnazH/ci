#define BOOST_TEST_MODULE "FCI_Davidson_test"


#include <hf.hpp>
#include <cpputil.hpp>

#include "FCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( FCI_H2_STO_3G_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation as reference
    ci::FCI fci_dense (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    fci_dense.solve(numopt::eigenproblem::SolverType::DENSE);


    // Do a Davidson FCI calculation
    ci::FCI fci_davidson (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    fci_davidson.solve(numopt::eigenproblem::SolverType::DAVIDSON);


    BOOST_CHECK(std::abs(fci_dense.get_eigenvalue() - fci_davidson.get_eigenvalue()) < 1.0e-12);
    BOOST_CHECK(cpputil::linalg::areEqualEigenvectors(fci_dense.get_eigenvector(), fci_davidson.get_eigenvector(), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( FCI_H2_6_31Gxx_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Prepare the AO basis
    libwint::Molecule h2 ("../tests/reference_data/h2_cristina.xyz");
    libwint::AOBasis ao_basis (h2, "6-31G**");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation as reference
    ci::FCI fci_dense (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    fci_dense.solve(numopt::eigenproblem::SolverType::DENSE);


    // Do a Davidson FCI calculation
    ci::FCI fci_davidson (so_basis, 1, 1);  // N_alpha = 1, N_beta = 1
    fci_davidson.solve(numopt::eigenproblem::SolverType::DAVIDSON);


    BOOST_CHECK(std::abs(fci_dense.get_eigenvalue() - fci_davidson.get_eigenvalue()) < 1.0e-12);
    BOOST_CHECK(cpputil::linalg::areEqualEigenvectors(fci_dense.get_eigenvector(), fci_davidson.get_eigenvector(), 1.0e-09));
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_STO_3G_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Prepare the AO basis
    libwint::Molecule h2o ("../tests/reference_data/h2o.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h2o, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a dense FCI calculation as reference
    ci::FCI fci_dense (so_basis, 5, 5);  // N_alpha = 5, N_beta = 5
    fci_dense.solve(numopt::eigenproblem::SolverType::DENSE);


    // Do a Davidson FCI calculation
    ci::FCI fci_davidson (so_basis, 5, 5);  // N_alpha = 5, N_beta = 5
    fci_davidson.solve(numopt::eigenproblem::SolverType::DAVIDSON);


    BOOST_CHECK(std::abs(fci_dense.get_eigenvalue() - fci_davidson.get_eigenvalue()) < 1.0e-12);
    BOOST_CHECK(cpputil::linalg::areEqualEigenvectors(fci_dense.get_eigenvector(), fci_davidson.get_eigenvector(), 1.0e-08));
}

BOOST_AUTO_TEST_CASE ( FCI_test ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Prepare the AO basis
    libwint::Molecule h4 ("../tests/h4_test.xyz");
    libwint::AOBasis ao_basis (h4, "6-31G**");
    ao_basis.calculateIntegrals();

    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (h4, ao_basis, 1.0e-06);
    rhf.solve();
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());

    // Do a Davidson FCI calculation
    ci::FCI fci_davidson (so_basis, 2, 2);
    fci_davidson.solve(numopt::eigenproblem::SolverType::DAVIDSON);

    std::cout << std::setprecision(15) << fci_davidson.get_eigenvalue() << std::endl;
}
