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
#include "FCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>



namespace ci {


/*
 *  OVERRIDDEN PRIVATE METHODS
 */

/**
 *  Given a @param matrix_solver, construct the FCI Hamiltonian matrix in the solver's matrix representation.
 */
void FCI::constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) {

    // The construction of the FCI Hamiltonian is implemented in three parts: alpha-alpha, beta-beta, and alpha-beta


    // 1. ALPHA-ALPHA
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_alpha (0, this->addressing_scheme_alpha);  // alpha spin string with address 0


    for (size_t I_alpha = 0; I_alpha < this->dim_alpha; I_alpha++) {  // I_alpha loops over all the addresses of the alpha spin strings

        size_t coupling_address_index = 0;  // index of |J_alpha> in the (N_alpha * (K + 1 - N_alpha))-long std::vector
                                            // located at alpha_one_electron_couplings[I_alpha]

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign for the annihilation operator (a_p)

            if (spin_string_alpha.annihilate(p, sign_p)) {
                for (size_t q = 0; q < this->K; q++) {  // q loops over SOs

                    // one-electron contributions for alpha, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_alpha.create(q, sign_pq)) {

                        size_t J_alpha = spin_string_alpha.address(this->addressing_scheme_alpha);

                        // For the 'diagonal beta contributions', i.e. I_beta = J_beta, the one-electron alpha contributions
                        // are the same
                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta
                        for (size_t I_beta = 0; I_beta < this->dim_beta; I_beta++) {
                            double value = sign_pq * this->so_basis.get_h_SO(p,q);
                            matrix_solver->addToMatrix(value, I_alpha * this->dim_beta + I_beta, J_alpha * this->dim_beta + I_beta);
                        }

                        // We have found a spin string that is one electron excitation away from |I_alpha>
                        // We will store it, since these strings are also needed in the alpha-beta part
                        this->alpha_one_electron_couplings[I_alpha][coupling_address_index] = OneElectronCoupling{sign_pq, p, q, J_alpha};
                        coupling_address_index++;
                        spin_string_alpha.annihilate(q);  // undo the previous creation on q
                    }  // create on q (alpha)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                                       // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_alpha.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < this->K; r++) {
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                                if (spin_string_alpha.create(r, sign_pqr)) {
                                for (size_t s = 0; s < this->K; s++) {

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_alpha.create(s, sign_pqrs)) {

                                        size_t Ja = spin_string_alpha.address(this->addressing_scheme_alpha);

                                        // For the 'diagonal beta contributions', i.e. Ib = Jb, the two-electron alpha
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_b + I_b
                                        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {
                                            double value = sign_pqrs * 0.5 * this->so_basis.get_g_SO(s,p,r,q);
                                            matrix_solver->addToMatrix(value, I_alpha * this->dim_beta + Ib, Ja * this->dim_beta + Ib);
                                        }

                                        spin_string_alpha.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (alpha)
                                }  // loop over s

                                spin_string_alpha.annihilate(r);  // undo the previous creation on r
                            }  // create on r (alpha)
                        }  // loop over r

                        spin_string_alpha.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (alpha)
                }  // loop over q

                spin_string_alpha.create(p);  // undo the previous annihilation on p
            }  // annihilate p (alpha)
        }  // loop over p


        if (I_alpha < this->dim_alpha - 1) {  // prevent the last permutation to occur
            spin_string_alpha.nextPermutation();
        }
    }  // loop over alpha addresses (I_alpha)


    // 2. BETA-BETA
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_beta (0, this->addressing_scheme_beta);  // beta spin string with address 0

    for (size_t I_beta = 0; I_beta < this->dim_beta; I_beta++) {  // I_beta loops over addresses of all beta spin strings

        size_t coupling_address_index = 0;  // index of |J_beta> in the (N_beta * (K + 1 - N_beta))-long std::vector
                                            // located at alpha_one_electron_couplings[I_alpha]

        for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;

            if (spin_string_beta.annihilate(p, sign_p)) {
                for (size_t q = 0; q < this->K; q++) {  // q loops over SOs

                    // one-electron contributions for beta, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_beta.create(q, sign_pq)) {

                        size_t J_beta = spin_string_beta.address(this->addressing_scheme_beta);

                        // For the 'diagonal alpha contributions', i.e. I_alpha = J_alpha, the one-electron beta contributions are
                        // the same
                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta

                        for (size_t I_alpha = 0; I_alpha < this->dim_alpha; I_alpha++) {
                            double value = sign_pq * this->so_basis.get_h_SO(p,q);
                            matrix_solver->addToMatrix(value, I_alpha * this->dim_beta + I_beta, I_alpha * this->dim_beta + J_beta);
                        }

                        // We have found a spin string that is one electron excitation away from |I_alpha>
                        // We will store it, since these strings are also needed in the alpha-beta part
                        this->beta_one_electron_couplings[I_beta][coupling_address_index] = OneElectronCoupling{sign_pq, p, q, J_beta};
                        coupling_address_index++;
                        spin_string_beta.annihilate(q);  // undo the previous creation on q
                    }  // create on q (beta)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                                       // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_beta.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < this->K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                            if (spin_string_beta.create(r, sign_pqr)) {
                                for (size_t s = 0; s < this->K; s++) {  // s loops over SOs

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_beta.create(s, sign_pqrs)) {

                                        size_t Jb = spin_string_beta.address(this->addressing_scheme_beta);

                                        // For the 'diagonal alpha contributions', i.e. Ia = Ja, the two-electron beta
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address IaIb = Ia * dim_b + I_b
                                        for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {
                                            double value = sign_pqrs * 0.5 * this->so_basis.get_g_SO(s,p,r,q);
                                            matrix_solver->addToMatrix(value, Ia * dim_beta + I_beta, Ia * dim_beta + Jb);
                                        }

                                        spin_string_beta.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (beta)
                                }  // loop over s

                                spin_string_beta.annihilate(r);  // undo the previous creation on r
                            }  // create on r (beta)
                        }  // loop over r

                        spin_string_beta.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (beta)
                }  // loop over q

                spin_string_beta.create(p);  // undo the previous annihilation on p
            } // annihilate on p (beta)
        }  // loop over p

        if (I_beta < dim_beta - 1) {  // prevent last permutation to occur
            spin_string_beta.nextPermutation();
        }
    }  // loop over beta addresses (I_beta)


    // 3. ALPHA-BETA
    for (size_t I_alpha = 0; I_alpha < this->dim_alpha; I_alpha++) {  // loop over alpha addresses
        for (size_t I_beta = 0; I_beta < this->dim_beta; I_beta++) {  // loop over beta addresses

            for (const auto& alpha : this->alpha_one_electron_couplings[I_alpha]) {  // traverse all OneElectronCouplings for I_alpha
                for (const auto& beta : this->beta_one_electron_couplings[I_beta]) {  // traverse all OneElectronCouplings for I_beta

                    int sign = alpha.sign * beta.sign;
                    double value = sign * this->so_basis.get_g_SO(alpha.p,alpha.q,beta.p,beta.q);
                    matrix_solver->addToMatrix(value, I_alpha * this->dim_beta + I_beta, alpha.address * this->dim_beta + beta.address);  // alpha is the major index
                }  // beta OneElectronCouplings
            }  // alpha OneElectronCouplings

        }  // loop over beta addresses (I_beta)
    }  // loop over alpha addresses (I_alpha)
}


/**
 *  @return the action of the FCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd FCI::matrixVectorProduct(const Eigen::VectorXd& x) {
    throw std::logic_error("This function hasn't been implemented yet.");
}


/**
 *  @set the diagonal of the matrix representation of the FCI Hamiltonian.
 */
void FCI::calculateDiagonal() {


    // Initialize the diagonal
    this->diagonal = Eigen::VectorXd::Zero(this->dim);


    // Calculate the effective one-electron integrals
    // TODO: move this to libwint
    Eigen::MatrixXd k_SO = this->so_basis.get_h_SO();
    for (size_t p = 0; p < this->K; p++) {
        for (size_t q = 0; q < this->K; q++) {
            for (size_t r = 0; r < this->K; r++) {
                k_SO(p,q) -= 0.5 * this->so_basis.get_g_SO(p,r,r,q);
            }
        }
    }


    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string_alpha (0, this->addressing_scheme_alpha);
    for (size_t Ia = 0; Ia < this->dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        bmqc::SpinString<unsigned long> spin_string_beta (0, this->addressing_scheme_beta);
        for (size_t Ib = 0; Ib < this->dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t p = 0; p < this->K; p++) {  // p loops over SOs

                if (spin_string_alpha.isOccupied(p)) {  // p is in Ia
                    this->diagonal(Ia * this->dim_beta + Ib) += k_SO(p, p);

                    for (size_t q = 0; q < this->K; q++) {  // q loops over SOs
                        if (spin_string_alpha.isOccupied(q)) {  // q is in Ia
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, p, q, q);
                        } else {  // q is not in I_alpha
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, q, q, p);
                        }

                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            this->diagonal(Ia * this->dim_beta + Ib) += this->so_basis.get_g_SO(p, p, q, q);
                        }
                    }  // q loop
                }


                if (spin_string_beta.isOccupied(p)) {  // p is in Ib
                    this->diagonal(Ia * this->dim_beta + Ib) += k_SO(p, p);


                    for (size_t q = 0; q < this->K; q++) {  // q loops over SOs
                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, p, q, q);

                        } else {  // q is not in I_beta
                            this->diagonal(Ia * this->dim_beta + Ib) += 0.5 * this->so_basis.get_g_SO(p, q, q, p);
                        }
                    }  // q loop
                }

            }  // p loop

            if (Ib < this->dim_beta - 1) {  // prevent last permutation to occur
                spin_string_beta.nextPermutation();
            }
        }  // beta address (Ib) loop

        if (Ia < this->dim_alpha - 1) {  // prevent last permutation to occur
            spin_string_alpha.nextPermutation();
        }
    }  // alpha address (Ia) loop
}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param so_basis, a number of alpha electrons @param N_alpha and a number of beta electrons
 *  @param N_beta.
 */
FCI::FCI(libwint::SOBasis& so_basis, size_t N_alpha, size_t N_beta) :
        BaseCI(so_basis, ci::FCI::calculateDimension(so_basis.get_K(), N_alpha, N_beta)),
        K (so_basis.get_K()),
        N_alpha (N_alpha),
        N_beta (N_beta),
        addressing_scheme_alpha (bmqc::AddressingScheme(this->K, this->N_alpha)),
        addressing_scheme_beta (bmqc::AddressingScheme(this->K, this->N_beta)),
        dim_alpha (ci::FCI::calculateDimension(so_basis.get_K(), N_alpha, 0)),
        dim_beta (ci::FCI::calculateDimension(so_basis.get_K(), 0, N_beta)),
        alpha_one_electron_couplings (this->dim_alpha, std::vector<OneElectronCoupling>(this->N_alpha * (this->K + 1 - this->N_alpha))),
        beta_one_electron_couplings (this->dim_beta, std::vector<OneElectronCoupling>(this->N_beta * (this->K + 1 - this->N_beta)))
{}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  Given a number of spatial orbitals @param K, a number of alpha electrons @param N_A, and a number of beta electrons
 *  @param N_B, @return the dimension of the FCI space.
 */
size_t FCI::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {

    // K, N_alpha, N_beta are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double_alpha = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_alpha));
    auto dim_double_beta = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_beta));
    auto dim_double = dim_double_alpha * dim_double_beta;

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  Calculate all the 1-RDMs.
 */
void FCI::calculate1RDMs() {
    //initialize as zero matrix.
    this->one_rdm_aa = Eigen::MatrixXd::Zero(this->K,this->K);
    this->one_rdm_bb = Eigen::MatrixXd::Zero(this->K,this->K);

    // Alpha-branch
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_alpha (0, this->addressing_scheme_alpha);  // spin string with address 0
    for (size_t I = 0; I < this->dim_alpha; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string_alpha.nextPermutation();
        }for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;
            if(spin_string_alpha.annihilate(p,sign_p)){
                double coefficient =  0;
                for(size_t i = 0;i<dim_beta;i++){  //repeat for each beta possibility
                    coefficient += this->eigensolver_ptr->get_eigenvector()(I*dim_beta+i)*this->eigensolver_ptr->get_eigenvector()(I*dim_beta+i);
                }  // dim_beta loop
                this->one_rdm_aa(p,p) += coefficient;
                for(size_t q = 0; q < p;q++){
                    int sign_q = sign_p;
                    if(spin_string_alpha.create(q,sign_q)){
                        size_t J = spin_string_alpha.address(this->addressing_scheme_alpha);

                        double coefficient = 0;
                        for(size_t i = 0;i<dim_beta;i++){
                            coefficient += this->eigensolver_ptr->get_eigenvector()(I*dim_beta+i)
                                           *this->eigensolver_ptr->get_eigenvector()(J*dim_beta+i);
                        }  // dim_beta loop
                        one_rdm_aa(q,p) += sign_q*coefficient;
                        one_rdm_aa(p,q) += sign_q*coefficient;

                        spin_string_alpha.annihilate(q);

                    }  // create q
                }  // q
                spin_string_alpha.create(p);
            }  // anni p
        }  // p
    }  // I

    // Beta -Branch
    bmqc::SpinString<boost::dynamic_bitset<>> spin_string_beta (0, this->addressing_scheme_beta);  // spin string with address 0
    for (size_t I = 0; I < this->dim_beta; I++) {  // I loops over all the addresses of the spin strings
        if (I > 0) {
            spin_string_beta.nextPermutation();
        }for (size_t p = 0; p < this->K; p++) {  // p loops over SOs
            int sign_p = 1;
            if(spin_string_beta.annihilate(p,sign_p)){

                double coefficient =  0;
                for(size_t i = 0;i<dim_alpha;i++){  //repeat for each alpha possibility
                    coefficient += this->eigensolver_ptr->get_eigenvector()(i*dim_beta+I)*this->eigensolver_ptr->get_eigenvector()(i*dim_beta+I);
                }  // dim_alpha loop
                this->one_rdm_bb(p,p) += coefficient;
                for(size_t q = 0; q < p;q++){
                    int sign_q = sign_p;
                    if(spin_string_beta.create(q,sign_q)){
                        size_t J = spin_string_beta.address(this->addressing_scheme_beta);

                        double coefficient = 0;
                        for(size_t i = 0;i<dim_alpha;i++){
                            coefficient += this->eigensolver_ptr->get_eigenvector()(i*dim_beta+I)
                                           *this->eigensolver_ptr->get_eigenvector()(i*dim_beta+J);
                        }
                        one_rdm_bb(q,p) += sign_q*coefficient;
                        one_rdm_bb(p,q) += sign_q*coefficient;
                        spin_string_beta.annihilate(q);

                    }  // create q


                }  // q
                spin_string_beta.create(p);
            }  // anni p
        }  // p
    }  // I

    this->one_rdm = this->one_rdm_aa + this->one_rdm_bb;
    this->are_computed_one_rdms = true;
}

/**
 *  Calculate all the 2-RDMs.
 */
void FCI::calculate2RDMs() {
    throw std::logic_error("This function hasn't been implemented yet.");

}



}  // namespace ci
