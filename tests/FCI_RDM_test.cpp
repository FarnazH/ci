#define BOOST_TEST_MODULE "FCI_RDM_test"



#include "FCI.hpp"

#include "RHF.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



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

/**
 *  @return the reduced-over 2-RDM @param d: D(p,q) = d(p,q,r,r)
 *
 *  In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
 *  Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves
 */
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



/*
 *  TESTS
 */
BOOST_AUTO_TEST_CASE ( H2O_1RDM_spin_trace_FCI ) {

    // Test if the traces of the spin-resolved 1-RDMs gives the number

    // Prepare an SOBasis
    size_t N_alpha = 5;
    size_t N_beta = 5;

    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Calculate the FCI 1-RDMs
    ci::FCI fci (so_basis, N_alpha, N_beta);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    fci.calculate1RDMs();


    Eigen::MatrixXd D_aa = fci.get_one_rdm_aa();
    BOOST_CHECK(std::abs(D_aa.trace() - N_alpha) < 1.0e-12);

    Eigen::MatrixXd D_bb = fci.get_one_rdm_bb();
    BOOST_CHECK(std::abs(D_bb.trace() - N_beta) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( H2O_2RDM_spin_trace_FCI ) {

    // Test if the traces of the spin-resolved 2-RDMs (d_ppqq) gives the correct number

    // Prepare the AO basis
    size_t N_alpha = 5;
    size_t N_beta = 5;

    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Do a dense FCI calculation based on a given SO basis, and calculate the 2-RDMs
    ci::FCI fci (so_basis, N_alpha, N_beta);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    fci.calculate2RDMs();


    BOOST_CHECK(std::abs(trace_2RDM(fci.get_two_rdm_aaaa()) - N_alpha*(N_alpha-1)) < 1.0e-12);
    BOOST_CHECK(std::abs(trace_2RDM(fci.get_two_rdm_aabb()) - N_alpha*N_beta) < 1.0e-12);
    BOOST_CHECK(std::abs(trace_2RDM(fci.get_two_rdm_bbaa()) - N_beta*N_beta) < 1.0e-12);
    BOOST_CHECK(std::abs(trace_2RDM(fci.get_two_rdm_bbbb()) - N_beta*(N_beta-1)) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( H2O_1RDM_2RDM_trace_FCI ) {

    // Test if the relevant 2-RDM trace reduction gives the 1-RDM for DOCI


    // Prepare the AO basis
    size_t N_alpha = 5;
    size_t N_beta = 5;
    size_t N = N_alpha + N_beta;

    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Do a dense FCI calculation based on a given SO basis, and calculate the 1- and 2-RDMs
    ci::FCI fci (so_basis, N_alpha, N_beta);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);
    fci.calculate1RDMs();
    fci.calculate2RDMs();


    Eigen::MatrixXd D = fci.get_one_rdm();
    Eigen::Tensor<double, 4> d = fci.get_two_rdm();

    Eigen::MatrixXd D_from_reduction = (1.0/(N-1)) * reduce_2RDM(d);

    BOOST_CHECK(D.isApprox(D_from_reduction, 1.0e-12));
}
