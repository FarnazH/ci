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
#ifndef CI_FCI_HPP
#define CI_FCI_HPP


#include "BaseCI.hpp"



namespace ci {



class FCI : public ci::BaseCI {
private:
    const size_t dim_alpha;  // the dimension of the alpha CI space
    const size_t dim_beta;  // the dimension of the beta CI space

    const size_t N_alpha;  // number of alpha electrons
    const size_t N_beta;  // number of beta electrons

    const bmqc::AddressingScheme addressing_scheme_alpha;
    const bmqc::AddressingScheme addressing_scheme_beta;

    // Rectangular matrix of SpinEvaluations
    /**
     *  A small struct that is used to hold in memory the @param address of spin strings differing in one electron
     *  excitation (an annihilation on orbital @param p and a creation on orbital @param q) that are coupled through the
     *  Hamiltonian.
     *
     *  During the construction of the FCI Hamiltonian, the one-electron excited coupling strings are both needed in the
     *  alpha, beta, and alpha-beta parts. When a spin string is found that couples to another spin string (with address
     *  I), the address of the coupling spin string is hold in memory, in the following way: in a
     *  std::vector<std::vector<OneElectronCoupling>> (with dimension I_alpha * N_alpha * (K + 1 - N_alpha)), at every outer index
     *  I_alpha, a std::vector of OneElectronCouplings is kept, each coupling through the Hamiltonian to that particular
     *  spin string with address I_alpha. Of course, the beta case is similar.
     *
     *  The @param sign of the matrix element, i.e. <I_alpha | H | address> is also stored as a parameter.
     *
     *
     *  We can keep this many addresses in memory because the resulting dimension (cfr. dim_alpha * N_alpha * (K + 1 - N_alpha)) is
     *  significantly less than the dimension of the FCI space (cfr. I_alpha * I_beta).
     *
     *  The number of coupling spin strings for an alpha string is equal to N_alpha * (K + 1 - N_alpha), since we have to pick
     *  one out of N_alpha occupied indices to annihilate, and afterwards (after the annihilation) we have (K + 1 - N_A)
     *  choices to pick an index to create on.
     */
    struct OneElectronCoupling {
        int sign;
        size_t p;
        size_t q;
        size_t address;
    };
    
    // The following are rectangular arrays of dimension (dim_alpha * N_alpha * (K + 1 - N_alpha)) and similarly for beta,
    // storing one-electron excited coupling addresses (cfr. the documentation about the OneElectronCoupling struct)
    std::vector<std::vector<OneElectronCoupling>> alpha_one_electron_couplings;  
    std::vector<std::vector<OneElectronCoupling>> beta_one_electron_couplings;



    // OVERRIDDEN PRIVATE METHODS
    /**
     *  Given a @param matrix_solver, construct the FCI Hamiltonian matrix in the solver's matrix representation.
     */
    void constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) override;

    /**
     *  @return the action of the FCI Hamiltonian of the coefficient vector @param x.
     */
    Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) override;

    /**
     *  @set the diagonal of the matrix representation of the FCI Hamiltonian.
     */
    void calculateDiagonal() override;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param so_basis, a number of alpha electrons @param N_alpha and a number of beta electrons
     *  @param N_beta.
     */
    FCI(libwint::SOBasis& so_basis, size_t N_alpha, size_t N_beta);


    // DESTRUCTOR
    ~FCI() override = default;


    // GETTERS
    size_t get_dim_alpha() const { return this->dim_alpha; };
    size_t get_dim_beta() const { return this->dim_beta; };


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K, a number of alpha electrons @param N_A, and a number of beta electrons
     *  @param N_B, @return the dimension of the FCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_alpha, size_t N_beta);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  Calculate all the 1-RDMs.
     */
    void calculate1RDMs() override;

    /**
     *  Calculate all the 2-RDMs.
     */
    void calculate2RDMs() override;
};


}  // namespace ci



#endif  // CI_FCI_HPP
