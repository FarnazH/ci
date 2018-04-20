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
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
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


BOOST_AUTO_TEST_CASE ( lih_1RDM_2RDM_contraction_DOCI ) {

    // Test if the relevant 2-RDM contraction gives the 1-RDM for DOCI


    // Get the 1- and 2-RDMs from DOCI
    size_t N = 4;  // 4 electrons
    libwint:: SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, N);
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    doci.calculate1RDMs();
    doci.calculate2RDMs();

    Eigen::MatrixXd D = doci.get_one_rdm();
    Eigen::Tensor<double, 4> d = doci.get_two_rdm();


    // Contract the 2-RDM
    //      Specify the dimension that should be 'reduced' over     1/(N-1) d(p q r r)
    Eigen::array<int, 2> reductions ({2,3});
    //      Perform the contraction
    Eigen::Tensor<double, 2> reduction = d.sum(reductions);
    std::cout <<  reduction << std::endl;


    std::cout << D << std::endl;


    // Test D(0,0)
    double D00 = 0.0;
    for (size_t r = 0; r < 16; r++) {
        D00 += d(0,0,r,r);
    }
    std::cout << D00 << std::endl;  //  this is correct by inspection with D(0,0)

    // Test D(0,1)
    double D01 = 0.0;
    for (size_t r = 0; r < 16; r++) {
        D01 += d(0,1,r,r);
    }
    std::cout << D01 << std::endl;  // this is correct by inspection with D(0,1)



}


BOOST_AUTO_TEST_CASE ( compute_test ) {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../tests/reference_data/beh_cation_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);
    doci.calculate1RDMs();
    doci.calculate2RDMs();

    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_CASE ( reference_tests ) {

    // TODO: add reference
    BOOST_CHECK(false);
}
