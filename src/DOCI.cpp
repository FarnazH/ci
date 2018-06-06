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
#include "DOCI.hpp"

#include <cpputil.hpp>

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>
#include <unsupported/Eigen/MatrixFunctions>


#include <chrono>  // I don't think we need this anymore



namespace ci {


/*
 *  OVERRIDDEN PRIVATE METHODS
 */

/**
 *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation.
 */
void DOCI::constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) {

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string(0, this->addressing_scheme);  // spin string with address 0



    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings

        // Diagonal contribution
        matrix_solver->addToMatrix(this->diagonal(I), I, I);

        // Off-diagonal contribution
        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.isOccupied(p)) {  // if p in I
                for (size_t q = 0; q < p; q++) {  // q loops over SOs
                    if (!spin_string.isOccupied(q)) {  // if q not in I

                        spin_string.annihilate(p);
                        spin_string.create(q);

                        size_t J = spin_string.address(this->addressing_scheme);  // J is the address of a string that couples to I

                        // The loops are p->K and q<p. So, we should normally multiply by a factor 2 (since the summand is symmetric)
                        // However, we are setting both of the symmetric indices of Hamiltonian, so no factor 2 is required
                        matrix_solver->addToMatrix(this->so_basis.get_g_SO(p, q, p, q), I, J);
                        matrix_solver->addToMatrix(this->so_basis.get_g_SO(p, q, p, q), J, I);

                        spin_string.annihilate(q);  // reset the spin string after previous creation
                        spin_string.create(p);  // reset the spin string after previous annihilation
                    }
                }  // q < p loop
            }
        }  // p loop

        spin_string.nextPermutation();
    }  // address (I) loop
}


/**
 *  @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd DOCI::matrixVectorProduct(const Eigen::VectorXd& x) {

//    auto start = std::chrono::high_resolution_clock::now();


    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);  // spin string with address 0


    // Diagonal contributions
    Eigen::VectorXd matvec = this->diagonal.cwiseProduct(x);


    // Off-diagonal contributions
    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings

        for (size_t p = 0; p < this->K; p++) {  // p loops over all SOs
            if (spin_string.isOccupied(p)) {  // p in I
                for (size_t q = 0; q < p; q++) {  // q loops over all SOs smaller than p
                    if (!spin_string.isOccupied(q)) {  // q not in I

                        spin_string.annihilate(p);
                        spin_string.create(q);

                        size_t J = spin_string.address(this->addressing_scheme);  // J is the address of a string that couples to I

                        matvec(I) += this->so_basis.get_g_SO(p,q,p,q) * x(J);  // off-diagonal contribution
                        matvec(J) += this->so_basis.get_g_SO(p,q,p,q) * x(I);  // off-diagonal contribution for q > p (not explicitly in sum)

                        spin_string.annihilate(q);  // reset the spin string after previous creation
                        spin_string.create(p);  // reset the spin string after previous annihilation
                    }
                } // q < p loop
            }
        }  // p loop

        spin_string.nextPermutation();
    }  // address (I) loop

//    auto stop = std::chrono::high_resolution_clock::now();
//
//    std::cout << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
//                      << " microseconds in matvec." << std::endl;

    return matvec;
}


/**
 *  @set the diagonal of the matrix representation of the DOCI Hamiltonian.
 */
void DOCI::calculateDiagonal() {


    // When calling this->solve twice for the Davidson solver, we need to produce consistent results: (re)set the diagonal to zero!
    this->diagonal = Eigen::VectorXd::Zero(this->dim);


    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);

    for (size_t I = 0; I < this->dim; I++) {  // I loops over addresses of spin strings

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.isOccupied(p)) {  // p is in I
                this->diagonal(I) += 2 * this->so_basis.get_h_SO(p,p) + this->so_basis.get_g_SO(p,p,p,p);

                for (size_t q = 0; q < p; q++) {  // q loops over SOs
                    if (spin_string.isOccupied(q)) {  // q is in I

                        // Since we are doing a restricted summation q<p, we should multiply by 2 since the summand argument is symmetric.
                        this->diagonal(I) += 2 * (2*this->so_basis.get_g_SO(p,p,q,q) - this->so_basis.get_g_SO(p,q,q,p));
                    }
                }  // q loop
            }
        }  // p loop

        spin_string.nextPermutation();
    }  // address (I) loop
}


/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis and a number of electrons @param N.
 */
DOCI::DOCI(libwint::SOBasis& so_basis, size_t N) :
    BaseCI(so_basis, this->calculateDimension(so_basis.get_K(), N / 2)),
    K (so_basis.get_K()),
    N_P (N / 2),
    addressing_scheme (bmqc::AddressingScheme(this->K, this->N_P))  // since in DOCI, alpha==beta, we should make an
                                                                    // addressing scheme with the number of PAIRS.
{
    // Do some input checks
    if ((N % 2) != 0) {
        throw std::invalid_argument("You gave an odd amount of electrons, which is not suitable for DOCI.");
    }

    // We can already calculate the diagonal, since this only has to be done once
    this->calculateDiagonal();
}


/**
 *  Constructor based on a given @param so_basis and a @param molecule.
 */
DOCI::DOCI(libwint::SOBasis& so_basis, const libwint::Molecule& molecule) :
    DOCI (so_basis, molecule.get_N())
{}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  Given a number of spatial orbitals @param K and a number of electron pairs @param N_P, @return the dimension of
 *  the DOCI space.
 */
size_t DOCI::calculateDimension(size_t K, size_t N_P) {

    // K and N_P are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_P));

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  Calculate all the 1-RDMS for DOCI.
 */
void DOCI::calculate1RDMs() {

    // The formulas for the DOCI 1-RDMs can be found in (https://github.com/lelemmen/electronic_structure)


    this->one_rdm_aa = Eigen::MatrixXd::Zero(this->K,this->K);

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long>spin_string (0, this->addressing_scheme);  // spin string with address 0


    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string.nextPermutation();
        }

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.isOccupied(p)) {  // if p is occupied in I
                double c_I = this->eigensolver_ptr->get_eigenvector(I);  // coefficient of the I-th basis vector
                this->one_rdm_aa(p,p) += std::pow(c_I, 2);
            }
        }
    }

    // For DOCI, we have an additional symmetry
    this->one_rdm_bb = this->one_rdm_aa;

    this->one_rdm = this->one_rdm_aa + this->one_rdm_bb;

    this->are_computed_one_rdms = true;
}


/**
 *  Calculate all the 2-RDMS for DOCI.
 */
void DOCI::calculate2RDMs(){

    // The formulas for the DOCI 2-RDMs can be found in (https://github.com/lelemmen/electronic_structure)


    this->two_rdm_aaaa = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->two_rdm_aaaa.setZero();
    this->two_rdm_aabb = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->two_rdm_aabb.setZero();

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);  // spin string with address 0


    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string.nextPermutation();
        }

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            if (spin_string.annihilate(p)) {  // if p is occupied in I

                double c_I = this->eigensolver_ptr->get_eigenvector(I);  // coefficient of the I-th basis vector
                double c_I_2 = std::pow(c_I, 2);  // square of c_I

                this->two_rdm_aabb(p,p,p,p) += c_I_2;

                for (size_t q = 0; q < p; q++) {  // q loops over SOs with an index smaller than p
                    if (spin_string.create(q)) {  // if q is not occupied in I
                        size_t J = spin_string.address(this->addressing_scheme);  // the address of the coupling string
                        double c_J = this->eigensolver_ptr->get_eigenvector(J);  // coefficient of the J-th basis vector

                        this->two_rdm_aabb(p,q,p,q) += c_I * c_J;
                        this->two_rdm_aabb(q,p,q,p) += c_I * c_J;  // since we're looping for q < p

                        spin_string.annihilate(q);  // reset the spin string after previous creation on q
                    }

                    else {  // if q is occupied in I
                        this->two_rdm_aaaa(p,p,q,q) += c_I_2;
                        this->two_rdm_aaaa(q,q,p,p) += c_I_2;  // since we're looping for q < p

                        this->two_rdm_aaaa(p,q,q,p) -= c_I_2;
                        this->two_rdm_aaaa(q,p,p,q) -= c_I_2;  // since we're looping for q < p

                        this->two_rdm_aabb(p,p,q,q) += c_I_2;
                        this->two_rdm_aabb(q,q,p,p) += c_I_2;  // since we're looping for q < p
                    }
                }
                spin_string.create(p);  // reset the spin string after previous annihilation on p
            }
        }
    }

    // For DOCI, we have additional symmetries
    this->two_rdm_bbbb = this->two_rdm_aaaa;
    this->two_rdm_bbaa = this->two_rdm_aabb;

    this->two_rdm = this->two_rdm_aaaa + this->two_rdm_aabb + this->two_rdm_bbaa + this->two_rdm_bbbb;

    this->are_computed_two_rdms = true;
}


/**
 *  Perform Newton-step-based orbital optimization
 */
void DOCI::orbitalOptimize(numopt::eigenproblem::DavidsonSolverOptions* davidson_solver_options_ptr) {


    // TODO: decide on wanting to implement orbital optimization for other solvers
    // TODO: decide on if we want to make a copy of the solver options: the user should know its solver options are getting overwritten
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options = *davidson_solver_options_ptr;
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;

    size_t orbital_optimization_iterations = 0;
    size_t number_of_maximum_orbital_optimization_iterations = 128;
    double convergence_threshold = 1.0e-08;


    // Solve the DOCI eigenvalue equation, using the options provided
//    this->solve(davidson_solver_options_ptr);
    this->solve(&dense_solver_options);
    double old_eigenvalue = this->get_lowest_eigenvalue();
    std::cout << "Original energy: " << old_eigenvalue << std::endl;


    bool is_converged = false;
    while (!is_converged) {


        // Calculate the 1- and 2-RDMs
        this->calculate1RDMs();
        this->calculate2RDMs();


        // Calculate the electronic gradient at kappa=0
        Eigen::MatrixXd F = this->so_basis.calculateGeneralizedFockMatrix(this->one_rdm, this->two_rdm);
        Eigen::MatrixXd gradient_matrix = 2 * (F - F.transpose());
        Eigen::Map<Eigen::VectorXd> gradient_vector (gradient_matrix.data(), gradient_matrix.cols()*gradient_matrix.rows());  // at kappa = 0
//        std::cout << "Gradient vector: " << std::endl << gradient_vector << std::endl << std::endl;


        // Calculate the electronic Hessian at kappa=0
        Eigen::Tensor<double, 4> W = this->so_basis.calculateSuperGeneralizedFockMatrix(this->one_rdm, this->two_rdm);
        Eigen::Tensor<double, 4> hessian_tensor (this->K, this->K, this->K, this->K);
        hessian_tensor.setZero();

        for (size_t p = 0; p < this->K; p++) {
            for (size_t q = 0; q < this->K; q++) {
                for (size_t r = 0; r < this->K; r++) {
                    for (size_t s = 0; s < this->K; s++) {
                        hessian_tensor(p,q,r,s) = W(p,q,r,s) - W(p,q,s,r) + W(q,p,s,r) - W(q,p,r,s) + W(r,s,p,q) - W(r,s,q,p) + W(s,r,q,p) - W(s,r,p,q);
                    }
                }
            }
        }
        Eigen::MatrixXd hessian_matrix = cpputil::linalg::toMatrix(hessian_tensor);  // at kappa = 0
        assert(hessian_matrix.isApprox(hessian_matrix.transpose()));

        Eigen::array<int, 4> shuffle {1, 0, 2, 3};
        Eigen::Tensor<double, 4> hessian_tensor_qp = hessian_tensor.shuffle(shuffle);
        Eigen::MatrixXd hessian_matrix_qp = cpputil::linalg::toMatrix(hessian_tensor_qp);

        assert(hessian_matrix.isApprox(-hessian_matrix_qp));


        // Try to calculate the Cholesky decomposition of the Hessian to check if it is positive definite
        Eigen::LLT<Eigen::MatrixXd> LLT_of_hessian (hessian_matrix);
        if (LLT_of_hessian.info() == Eigen::NumericalIssue)  // the Hessian is not positive definite
        {

            // If the Hessian is not positive definite, we can try to modify it
            // Try to add a multiple (tau) of the identity: tau = - lowest eigenvalue of the Hessian
            // If we can't diagonalize, i.e. for large systems, we can use an algorithm to keep increasing tau until the Hessian is positive definite
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_solver (hessian_matrix);
            double lowest_eigenvalue = hessian_solver.eigenvalues()(0);
            std::cout << "lowest eigenvalue: " << lowest_eigenvalue << std::endl;

            hessian_matrix += -1.1 * lowest_eigenvalue * Eigen::MatrixXd::Identity(this->K*this->K, this->K*this->K);

            Eigen::LLT<Eigen::MatrixXd> LLT_of_modified_hessian (hessian_matrix);
            if (LLT_of_modified_hessian.info() == Eigen::NumericalIssue) {
                throw std::runtime_error("The modified Hessian still isn't invertible.");
            }
        }


        // Perform a Newton-step to find orbital rotation parameters kappa
        numopt::GradientFunction gradient_function = [gradient_vector](const Eigen::VectorXd& x) { return gradient_vector; };
        numopt::JacobianFunction hessian_function = [hessian_matrix](const Eigen::VectorXd& x) { return hessian_matrix; };


        //        Eigen::VectorXd kappa_vector = numopt::newtonStep(Eigen::VectorXd::Zero(this->K), gradient_function, hessian_function);
        Eigen::VectorXd kappa_vector = - hessian_matrix.inverse() * gradient_vector;
        // Change kappa back to a matrix
        Eigen::Map<Eigen::MatrixXd> kappa_matrix (kappa_vector.data(), this->K, this->K);


        if (!(hessian_matrix * kappa_vector + gradient_vector).isZero()) {  // the Hessian was singular, so use gradient descent
            std::cout << "Hessian wasn't invertible... Had to do gradient descent ..." << std::endl;
            kappa_matrix = -gradient_matrix;
        }

//        assert((kappa_vector_newton + hessian_matrix.inverse() * gradient_vector).isZero());


//        std::cout << "Kappa vector: " << std::endl << kappa_vector << std::endl << std::endl;


        // What if the Hessian is singular? Do a gradient descent  // TODO put this step into numopt
//        Eigen::MatrixXd kappa_matrix_descent = -gradient_matrix;


        if (kappa_matrix.norm() < convergence_threshold) {
            is_converged = true;
        } else {
            orbital_optimization_iterations++;

            // TODO: should I do something else?
            // The current solution can now be accessed using the normal getters

            if (orbital_optimization_iterations >= number_of_maximum_orbital_optimization_iterations ) {
                throw std::runtime_error("The OO-DOCI procedure failed to converge in the maximum number of allowed iterations.");
            }
        }


        std::cout << "Kappa matrix: " << std::endl << kappa_matrix << std::endl << std::endl;
        // Calculate the unitary rotation matrix
        Eigen::MatrixXd U = (kappa_matrix).exp();
        std::cout << "Unitary transformation matrix: " << std::endl << U << std::endl << std::endl;


        // Transform the integrals to the new orthonormal basis
        this->so_basis.rotate(U);


        //
        Eigen::VectorXd v = this->get_lowest_eigenvector();
//        std::cout << "Eigenvector of the old Hamiltonian: v: " << std::endl << v << std::endl << std::endl;
        Eigen::VectorXd Hv = this->matrixVectorProduct(v);
//        std::cout << "New Hamiltonian action on the old eigenvector: Hv: " << std::endl << Hv << std::endl << std::endl;

        std::cout << "Energy of the old eigenvector in the new Hamiltonian: " << std::endl << v.dot(Hv) << std::endl << std::endl;


        // Solve the DOCI eigenvalue problem in the new basis
        davidson_solver_options.X_0 = v;
//        this->solve(davidson_solver_options_ptr);
        this->solve(&dense_solver_options);
    }  // while not converged
}


}  // namespace ci
