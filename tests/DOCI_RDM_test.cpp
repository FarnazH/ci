#define BOOST_TEST_MODULE "DOCI_RDM_test"



#include "DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( lih_energy_RDM_contraction_DOCI ) {

    // Test if the contraction of the 1- and 2-RDMs with the one- and two-electron integrals gives the DOCI energy


    // Get the DOCI energy as the lowest eigenvalue of the dense DOCI Hamiltonian
    libwint:: SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    Eigen::MatrixXd h = so_basis.get_h_SO();
    Eigen::Tensor<double, 4> g = so_basis.get_g_SO();

    ci::DOCI doci (so_basis, 4);  // 4 electrons
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci.solve(&dense_options);
    double energy_by_eigenvalue = doci.get_eigenvalue();


    // Calculate the DOCI energy as the relevant contraction with the one- and two-electron integrals
    doci.calculate1RDMs();
    Eigen::MatrixXd D = doci.get_one_rdm();
    double energy_by_contraction = (h * D).trace();

    doci.calculate2RDMs();
    Eigen::Tensor<double, 4> d = doci.get_two_rdm();

    // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
    //      0.5 g(p q r s) d(p q r s)
    Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
    //      Perform the contraction
    Eigen::Tensor<double, 0> contraction = 0.5 * g.contract(d, contractions);

    // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
    energy_by_contraction += contraction(0);


    BOOST_CHECK(std::abs(energy_by_eigenvalue - energy_by_contraction) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_1RDM_2RDM_trace_DOCI ) {

    // Test if the relevant 2-RDM trace gives the 1-RDM for DOCI


    // Get the 1- and 2-RDMs from DOCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    libwint:: SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", K);
    ci::DOCI doci (so_basis, N);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci.solve(&dense_options);
    doci.calculate1RDMs();
    doci.calculate2RDMs();

    Eigen::MatrixXd D = doci.get_one_rdm();
    Eigen::Tensor<double, 4> d = doci.get_two_rdm();


    // Trace the 2-RDM
    //      Specify the dimension that should be 'reduced' over     d(p q r r)
    //      In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction
    // Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves
    Eigen::MatrixXd D_from_reduction = Eigen::MatrixXd::Zero(K, K);
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                D_from_reduction(p,q) += (1.0/(N-1)) * d(p,q,r,r);
            }
        }
    }

    BOOST_CHECK(D.isApprox(D_from_reduction, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( lih_1RDM_trace ) {

    // Test if the trace of the 1-RDM gives N


    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    libwint:: SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", K);
    ci::DOCI doci (so_basis, N);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci.solve(&dense_options);
    doci.calculate1RDMs();

    Eigen::MatrixXd D = doci.get_one_rdm();

    BOOST_CHECK(std::abs(D.trace() - N) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_2RDM_trace ) {

    // Test if the trace of the 2-RDM (d_ppqq) gives N(N-1)


    // Get the 2-RDM from DOCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    libwint:: SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", K);
    ci::DOCI doci (so_basis, N);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    doci.solve(&dense_options);
    doci.calculate2RDMs();

    Eigen::Tensor<double, 4> d = doci.get_two_rdm();

    // Trace the 2-RDM
    //      Specify the dimension that should be 'reduced' over     d(p p q q)
    //      In this case, the 'trace' tensor operation should be used, but the current Eigen3 (3.3.4) hasn't released that support yet
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction
    // Since Eigen3 hasn't released tensor.trace() yet, we will do the reduction ourselves
    double trace_value = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace_value += d(p,p,q,q);
        }
    }

    BOOST_CHECK(std::abs(trace_value - N*(N-1)) < 1.0e-12);
}
