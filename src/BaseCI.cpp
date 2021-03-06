// This file is part of GQCG-ci.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-ci is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-ci is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-ci.  If not, see <http://www.gnu.org/licenses/>.
#include "BaseCI.hpp"

#include <numopt.hpp>

#include "DenseSolver.hpp"
#include "SparseSolver.hpp"


#include <iomanip>


namespace ci {



/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  Protected constructor given a @param so_basis and a dimension @dim.
 */
BaseCI::BaseCI(libwint::SOBasis& so_basis, size_t dim) :
    so_basis (so_basis),
    K (this->so_basis.get_K()),
    dim (dim),
    diagonal (Eigen::VectorXd::Zero(this->dim))
{}



/*
 *  PROTECTED METHODS
 */

/**
 *  Initialize and subsequently solve the eigenvalue problem associated to the derived CI class, using a pointer
 *  to a @param matrix_solver_ptr.
 */
void BaseCI::solveMatrixEigenvalueProblem(numopt::eigenproblem::BaseMatrixSolver* matrix_solver_ptr) {

    // Initialize the Hamiltonian matrix and solve the eigenvalue problem associated to it.
    this->constructHamiltonian(matrix_solver_ptr);
    matrix_solver_ptr->solve();
}



/*
 *  DESTRUCTOR
 */

BaseCI::~BaseCI() {
    delete this->eigensolver_ptr;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Providing a @param solver_options_ptr, find (the) lowest energy eigenpair(s) of the Hamiltonian
 */
void BaseCI::solve(numopt::eigenproblem::BaseSolverOptions* solver_options_ptr) {

    // The diagonal should be recalculated every time solve is called, since this->so_basis can have changed in between
    // subsequent diagonalization calls
    this->calculateDiagonal();


    // Depending on how the user wants to solve the eigenvalue problem, construct the appropriate solver
    switch (solver_options_ptr->get_solver_type()) {

        case numopt::eigenproblem::SolverType::DENSE: {
            auto dense_solver = new numopt::eigenproblem::DenseSolver(this->dim);  // deleted in the destructor
            this->solveMatrixEigenvalueProblem(dense_solver);
            this->eigensolver_ptr = dense_solver;  // prevent data from going out of scope
                                                   // we are only assigning this->eigensolver_ptr now, because
                                                   // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

        case numopt::eigenproblem::SolverType::SPARSE: {
            auto sparse_solver = new numopt::eigenproblem::SparseSolver(this->dim);  // deleted in the destructor
            this->solveMatrixEigenvalueProblem(sparse_solver);
            this->eigensolver_ptr = sparse_solver;  // prevent data from going out of scope
                                                    // we are only assigning this->eigensolver_ptr now, because
                                                    // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

        case numopt::eigenproblem::SolverType::DAVIDSON: {

            numopt::VectorFunction matrixVectorProduct = [this](const Eigen::VectorXd& x) { return this->matrixVectorProduct(x); };

            // Dynamic-cast to the derived pointer type, to be able to access X_0 (the initial guess(es))
            auto davidson_solver_options_ptr = dynamic_cast<numopt::eigenproblem::DavidsonSolverOptions*> (solver_options_ptr);  // it's OK if this pointer goes out of scope, we still have solver_options_ptr that's pointing to it

            this->eigensolver_ptr = new numopt::eigenproblem::DavidsonSolver(matrixVectorProduct, this->diagonal, davidson_solver_options_ptr->X_0);  // deleted in the destructor
            this->eigensolver_ptr->solve();
            break;
        }
    }
}



/*
 * GETTERS
 */

Eigen::MatrixXd BaseCI::get_one_rdm_aa() const {
    if (!this->are_computed_one_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this->one_rdm_aa;
}


Eigen::MatrixXd BaseCI::get_one_rdm_bb() const {
    if (!this->are_computed_one_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this-> one_rdm_bb;
}


Eigen::MatrixXd BaseCI::get_one_rdm() const {
    if (!this->are_computed_one_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this-> one_rdm;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_aaaa() const {
    if (!this->are_computed_two_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this->two_rdm_aaaa;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_aabb() const {
    if (!this->are_computed_two_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this->two_rdm_aabb;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_bbaa() const {
    if (!this->are_computed_two_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this->two_rdm_bbaa;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm_bbbb() const {
    if (!this->are_computed_two_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this->two_rdm_bbbb;
}


Eigen::Tensor<double, 4> BaseCI::get_two_rdm() const {
    if (!this->are_computed_two_rdms) {
        throw std::logic_error("The requested reduced density matrix is not computed yet.");
    }
    return this->two_rdm;
}



}  // namespace ci
