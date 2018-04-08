#include "DOCI.hpp"



#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>



#include <chrono>


void bmqc2_nextPermutation(unsigned long& spin_string) {

    // t gets this->representation's least significant 0 bits set to 1
    unsigned long t = spin_string | (spin_string - 1UL);

    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    spin_string = (t + 1UL) | (((~t & (t+1UL)) - 1UL) >> (__builtin_ctzl(spin_string) + 1UL));
}


size_t get_address(const unsigned long& spin_string, const bmqc::AddressingScheme& addressing_scheme, size_t start, size_t hits) {

    size_t copy = spin_string >> start;

    // An implementation of the formula in Helgaker, starting the addressing count from zero
    size_t address = 0;
    while(copy !=0 ){
        unsigned long t = copy & -copy;
        int r = __builtin_ctzl(copy);
        hits++;
        address += addressing_scheme.get_vertex_weights(start+r,hits);
        copy ^= t;
    }
    return address;
}


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

    auto start = std::chrono::high_resolution_clock::now();


    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    // TODO: determine when to switch from unsigned to unsigned long, unsigned long long or boost::dynamic_bitset<>
    bmqc::SpinString<unsigned long> spin_string_pre (0, this->addressing_scheme);  // spin string with address 0
    unsigned long spin_string = spin_string_pre.get_representation();

    // Diagonal contributions
    Eigen::VectorXd matvec = this->diagonal.cwiseProduct(x);


    // Off-diagonal contributions
    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        // We are essentially looking for all addresses J of spin strings that couple to the spin string at I
        // These are the ones in which:
        //      - an electron in an occupied orbital p is annihilated
        //      - an electron in an unoccupied orbital q is created

        double I_matvec_value = 0;  // accumulating in an auxiliary variable saves us the extra time associated to "matvec(I) +=" calls



        // We are dividing the address calculation in three parts:
        //      - address1:     the part of the address of the coupling spin string |J> with indices < p
        //      - address2:
        //      - address3:
        size_t counter1 = 0;  // counts the electrons in orbitals with indices < p in the coupling spin string |J>
        size_t address1 = 0;  // the part of the address of the coupling spin string |J> with indices < p


        // The first loop essentially iterates over all set bits in the spin string (i.e. indices p we can annihilate on)
        // We're going to keep in-place annihilating the right-most set bit, until we end up with no more set bits (i.e. "0")
        unsigned long copy = spin_string;
        while (copy != 0) {
            size_t p = __builtin_ctzl(copy);  // # of trailing zeros = index of the annihilation operator (i.e. must apply on index p)

            size_t address2 = 0;  // the part of the address of the spin string corresponding to the "gap"
            size_t gap = p+1; // Current position in the addressing scheme (called gap because
            size_t counter2 = counter1;  // counting electrons in the second loop)


            // Iterate over all unset bits that are beyond p (allowing us to calculate only the upper diagonal contributions (i.e. q > p))
            // If we want to use __builtin_ctz(), we must invert all bits with index > q
            unsigned long inverted = ~(copy | (copy-1));  // propagate the right-most set bit and invert the result
            while (__builtin_ctzl(inverted) < this->K) {  // we can't use inverted != 0 since the inverting also affects the bits > K
                size_t q = __builtin_ctzl(inverted);  // # of trailing zeros in the inverted spin string = index of the creation operator

                while (gap < q) { // allows us to bridge the gap (0's that were previously set bits) between two creation operators
                    // This allows us to update the address accordingly.
                    counter2++;
                    address2 += addressing_scheme.get_vertex_weights(gap,counter2);
                    gap++;
                }
                gap++; // set gap to the current position in the addressing scheme.


                size_t address3 = addressing_scheme.get_vertex_weights(q,counter2+1) + get_address(spin_string,addressing_scheme,q,counter2+1);  // the part of the address of the spin string with indices >= q


                // The final address J is is the sum of the three parts address1, address2 and address3
                size_t J = address1 + address2 + address3;

                I_matvec_value += this->so_basis.get_g_SO(p,q,p,q) * x(J);  // accumulating in an auxiliary variable saves us the extra time associated to matvec(I)+= calls
                matvec(J) += this->so_basis.get_g_SO(p,q,p,q) * x(I);

                inverted ^= inverted & -inverted;  // annihilate the least significant bit, i.e. create on the least significant position q>p in the original spin string
                                                   // the next unset bit in the original spin string can then again be found using __builtin_ctz(inverted)
            }  // q loop (creation)


            // Once all possible excitations (i.e. the suitable orbitals q>p) from an orbital p are found, we continue
            // with the next occupied orbital. Since we are destroying our "copy", we have to keep track of how many
            // electrons we have already encountered, since these contribute to the total address.
            counter1++;
            address1 += addressing_scheme.get_vertex_weights(p,counter1);
            copy ^= (copy & -copy);  // annihilate the least significant bit, i.e. on index p
                                     // the next set bit can then be again found using __builtin_ctz(copy)
        }  // p loop (annihilation)

        matvec(I) += I_matvec_value;
        bmqc2_nextPermutation(spin_string);
    }  // address (I) loop



    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                      << " microseconds in matvec." << std::endl;

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


}  // namespace ci
