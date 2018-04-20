#include "DOCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>



#include <chrono>



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
    addressing_scheme (bmqc::AddressingScheme(this->K, this->N_P)) // since in DOCI, alpha==beta, we should make an
                                                                   // addressing scheme with the number of PAIRS.
{
    // Do some input checks
    if ((N % 2) != 0) {
        throw std::invalid_argument("You gave an odd amount of electrons, which is not suitable for DOCI.");
    }

    if (this->K < this->N_P) {
        throw std::invalid_argument("Too many electrons to place into the given number of spatial orbitals.");
    }
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
    bmqc::SpinString<unsigned long> spin_string (0, this->addressing_scheme);  // spin string with address 0


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



}  // namespace ci
