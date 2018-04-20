#ifndef CI_DOCI_HPP
#define CI_DOCI_HPP



#include "BaseCI.hpp"



namespace ci {


class DOCI : public ci::BaseCI {
private:
    const size_t K;  // number of spatial orbitals
    const size_t N_P;  // number of electron pairs
    const bmqc::AddressingScheme addressing_scheme;


    // OVERRIDDEN PRIVATE METHODS
    /**
     *  Given a @param matrix_solver, construct the DOCI Hamiltonian matrix in the solver's matrix representation.
     */
    void constructHamiltonian(numopt::eigenproblem::BaseMatrixSolver* matrix_solver) override;

    /**
     *  @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
     */
    Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) override;

    /**
     *  @set the diagonal of the matrix representation of the Hamiltonian.
     */
    void calculateDiagonal() override;



public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param so_basis and a number of electrons @param N.
     */
    DOCI(libwint::SOBasis& so_basis, size_t N);

    /**
     *  Constructor based on a given @param so_basis and a @param molecule.
     */
    DOCI(libwint::SOBasis& so_basis, const libwint::Molecule& molecule);


    // DESTRUCTOR
    ~DOCI() override = default;


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K and a number of electron pairs @param N_P, @return the dimension of
     *  the DOCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_P);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  Calculate all the 1-RDMs for DOCI.
     */
    void calculate1RDMs() override;

    /**
     *  Calculate all the 2-RDMS for DOCI.
     */
    void calculate2RDMs() override;
};


}  // namespace ci



#endif  // CI_DOCI_HPP
