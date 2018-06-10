#define BOOST_TEST_MODULE "FCI_RDM_test"



#include "FCI.hpp"

#include "RHF.hpp"

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


BOOST_AUTO_TEST_CASE ( h2o_2RDM_aaaa_trace_fci ) {

    // Test if the trace of the aaaa-2-RDM (d_ppqq) gives N_alpha(N_alpha-1)

    // Prepare the AO basis
    size_t N_alpha = 5;
    size_t N_beta = 5;
    size_t K = 7;
    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, N_alpha, N_beta);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    fci.calculate2RDMs();


    Eigen::Tensor<double, 4> d_aaaa = fci.get_two_rdm_aaaa();

    // Trace the 2-RDM
    //      Specify the dimension that should be 'reduced' over     d(p p q q)
    //      In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction
    // Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves
    double trace_value = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace_value += d_aaaa(p,p,q,q);
        }
    }


    BOOST_CHECK(std::abs(trace_value - N_alpha*(N_alpha-1)) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( h2o_2RDM_bbbb_trace_fci ) {

    // Test if the trace of the bbbb-2-RDM (d_ppqq) gives N_beta(N_beta-1)

    // Prepare the AO basis
    size_t N_alpha = 5;
    size_t N_beta = 5;
    size_t K = 7;
    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, N_alpha, N_beta);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    fci.calculate2RDMs();


    Eigen::Tensor<double, 4> d_bbbb = fci.get_two_rdm_bbbb();

    // Trace the 2-RDM
    //      Specify the dimension that should be 'reduced' over     d(p p q q)
    //      In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction
    // Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves
    double trace_value = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace_value += d_bbbb(p,p,q,q);
        }
    }


    BOOST_CHECK(std::abs(trace_value - N_beta*(N_beta-1)) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( h2o_2RDM_aabb_trace_fci ) {

    // Test if the trace of the aabb-2-RDM (d_ppqq) gives N_alpha*N_beta

    // Prepare the AO basis
    size_t N_alpha = 5;
    size_t N_beta = 5;
    size_t K = 7;
    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, N_alpha, N_beta);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    fci.calculate2RDMs();


    Eigen::Tensor<double, 4> d_aabb = fci.get_two_rdm_aabb();

    // Trace the 2-RDM
    //      Specify the dimension that should be 'reduced' over     d(p p q q)
    //      In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction
    // Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves
    double trace_value = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace_value += d_aabb(p,p,q,q);
        }
    }


    BOOST_CHECK(std::abs(trace_value - N_alpha*N_beta) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( h2o_2RDM_bbaa_trace_fci ) {

    // Test if the trace of the bbaa-2-RDM (d_ppqq) gives N_alpha*N_beta

    // Prepare the AO basis
    size_t N_alpha = 5;
    size_t N_beta = 5;
    size_t K = 7;
    libwint::Molecule h2o ("../tests/reference_data/h2o_Psi4_GAMESS.xyz");
    libwint::AOBasis ao_basis (h2o, "STO-3G");
    ao_basis.calculateIntegrals();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Do a dense FCI calculation based on a given SO basis
    ci::FCI fci (so_basis, N_alpha, N_beta);
    fci.solve(numopt::eigenproblem::SolverType::DENSE);

    fci.calculate2RDMs();


    Eigen::Tensor<double, 4> d_bbaa = fci.get_two_rdm_bbaa();

    // Trace the 2-RDM
    //      Specify the dimension that should be 'reduced' over     d(p p q q)
    //      In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction
    // Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves
    double trace_value = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace_value += d_bbaa(p,p,q,q);
        }
    }


    BOOST_CHECK(std::abs(trace_value - N_beta*N_alpha) < 1.0e-12);
}
