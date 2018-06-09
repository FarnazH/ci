#define BOOST_TEST_MODULE "FCI_RDM_test"



#include "FCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( h2o_1RDM_trace_fci ) {

    // Test if the trace of the 1-RDM gives the number of electrons N

    // Prepare an SOBasis by transforming a non-orthogonal AOBasis with the Löwdin orthogonalization
    size_t N = 10;  // 10 electrons
    size_t K = 7;  // 7 SOs
    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    ci::FCI fci (so_basis, N/2, N-N/2);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    fci.calculate1RDMs();

    Eigen::MatrixXd D = fci.get_one_rdm();


    BOOST_CHECK(std::abs(D.trace() - N) < 1.0e-12);
}
