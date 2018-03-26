/**
 *  Do some timings concerning for dense and sparse DOCI calculations
 */


#include <chrono>
#include <iostream>

#include <Eigen/Dense>
#include <cpputil.hpp>

#include <numopt.hpp>
#include <libwint.hpp>
#include <ci.hpp>



/*
 *  HELPER FUNCTIONS
 */


/**
 *  Solve the dense DOCI problem (full diagonalization of the DOCI Hamiltonian matrix) for LiH.
 */
void solveDenseLiH() {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DENSE);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
}


/**
 *  Solve the sparse DOCI problem (find the lowest eigenpair of the DOCI Hamiltonian matrix) for LiH.
 */
void solveSparseLiH() {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::SPARSE);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
}


/**
 *  Solve the DOCI problem for LiH using Davidson's algorithm.
 */
void solveDavidsonLiH() {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
    ci::DOCI doci (so_basis, 4);  // 4 electrons
    doci.solve(numopt::eigenproblem::SolverType::DAVIDSON);


    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
}


/**
 *  Solve the DOCI problem for CO using Davidson's algorithm.
 */
void solveDavidsonCO() {

    // Do a DOCI calculation based on a given FCIDUMP file
    libwint::SOBasis so_basis ("../../tests/reference_data/co_631g_klaas.FCIDUMP", 28);  // 28 SOs
    ci::DOCI doci (so_basis, 14);  // 14 electrons
    doci.solve(numopt::eigenproblem::SolverType::DAVIDSON);


    // Calculate the total energy
    double internuclear_repulsion_energy = 2.2141305786610879e+01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
}


/**
 *  Print how long it takes to solve the dense DOCI problem for LiH.
 */
void printDenseLiHTimings() {

    cpputil::printExecutionTime("LiH dense DOCI", solveDenseLiH);  // (void *)() is implicitly converted to std::function<void ()>
}


/**
 *  Print how long it takes to solve the sparse DOCI problem for LiH.
 */
void printSparseLiHTimings() {

    cpputil::printExecutionTime("LiH sparse DOCI", solveSparseLiH);  // (void *)() is implicitly converted to std::function<void ()>
}


/**
 *  Print how long it takes to solve the DOCI problem for LiH with Davidson's algorithm.
 */
void printDavidsonLiHTimings() {

    cpputil::printExecutionTime("LiH Davidson DOCI", solveDavidsonLiH);  // (void *)() is implicitly converted to std::function<void ()>
}


/**
 *  Print how long it takes to solve the DOCI problem for CO with Davidson's algorithm.
 */
void printDavidsonCOTimings() {

    cpputil::printExecutionTime("CO Davidson DOCI", solveDavidsonCO);  // (void *)() is implicitly converted to std::function<void ()>
}


/*
 *  MAIN FUNCTION
 */

int main () {

//    printDenseLiHTimings();
//    printSparseLiHTimings();
//    printDavidsonLiHTimings();

    printDavidsonCOTimings();
}
