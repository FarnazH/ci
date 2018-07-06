#ifndef CI_ONVEXPANSION_HPP
#define CI_ONVEXPANSION_HPP


#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/range/adaptors.hpp>

#include "ONVExpansionComponent.hpp"
#include "FCI.hpp"
#include "DOCI.hpp"

#include <cpputil.hpp>



namespace ci {



/**
 *  A class that holds a wave function expansion: it is a wrapper around a std::vector<ci::ONVExpansionComponent<T>>,
 *  which holds alpha-ONVs, beta-ONVs and a corresponding coefficient
 *
 *  // TODO: should we make this unordered?
 *  // TODO: should we normalize / require normalization?
 */
template <typename T>
class ONVExpansion {
private:
    std::vector<ci::ONVExpansionComponent<T>> expansion;
    size_t K;  // number of spatial orbitals
    size_t dim;  // the dimension of the expansion


    bool are_calculated_one_rdms = false;
    bool are_calculated_two_rdms = false;

    Eigen::MatrixXd one_rdm_aa;  // alpha-alpha (a-a) 1-RDM
    Eigen::MatrixXd one_rdm_bb;  // beta-beta (b-b) 1-RDM
    Eigen::MatrixXd one_rdm;  // spin-summed (total) 1-RDM

    Eigen::Tensor<double, 4> two_rdm_aaaa;  // a-a-a-a 2-RDM
    Eigen::Tensor<double, 4> two_rdm_aabb;  // a-a-b-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm_bbaa;  // b-a-a-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm_bbbb;  // b-b-b-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm;  // spin-summed (total) 2-RDM

public:

    /*
     *  CONSTRUCTORS
     */
    /**
     *  Constructor based on a given initializer list
     */
    ONVExpansion(const std::initializer_list<ci::ONVExpansionComponent<T>>& list) :
        expansion (list),
        K (this->expansion[0].alpha.get_K()),
        dim (this->expansion.size())
    {
        // Check if all the ONVs have the same number of orbitals
        for (const auto& expansion_component : this->expansion) {
            if ((expansion_component.alpha.get_K() != this->K) || (expansion_component.beta.get_K() != this->K)) {
                throw std::invalid_argument("ONVExpansion(): The number of orbitals in the expansion is not consistent.");
            }
        }
        // TODO: Check for no duplicates?
    }


    /**
     *  Constructor based on a given @param GAMESS_filename
     */
    explicit ONVExpansion(const std::string& GAMESS_filename) {

        // If the filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
        std::ifstream input_file_stream(GAMESS_filename);

        if (!input_file_stream.good()) {
            throw std::runtime_error("ONVExpansion(): The provided GAMESS file is illegible. Maybe you specified a wrong path?");
        }


        // Do the actual parsing

        // Read in dummy lines up until we actually get to the ONVs and coefficients
        std::string line;
        while (std::getline(input_file_stream, line)) {
            if (line.find("ALPHA") != std::string::npos
                && line.find("BETA") != std::string::npos
                && line.find("COEFFICIENT") != std::string::npos) {  // if find returns an index that's different from the 'not-found' index

                std::getline(input_file_stream, line);  // this line should have dashes

                // We can actually get the number of orbitals by counting the number of dashes
                std::vector<std::string> splitted_line;
                boost::split(splitted_line, line, boost::is_any_of("|"));  // split on '|'

                this->K = splitted_line[0].length();
                break;
            }
        }


        // Read in the ONVs and the coefficients by splitting the line on '|', and then trimming whitespace
        // In the GAMESS format, the bit strings are ordered in reverse
        while (std::getline(input_file_stream, line)) {
            std::vector<std::string> splitted_line;
            boost::split(splitted_line, line, boost::is_any_of("|"));  // split on '|'


            // Create an alpha SpinString for the first field
            std::string trimmed_alpha = boost::algorithm::trim_copy(splitted_line[0]);
            if (trimmed_alpha.length() != this->K) {
                throw std::invalid_argument("ONVExpansion(): One of the provided alpha ONVs does not have the correct number of orbitals.");
            }
            std::string reversed_alpha (trimmed_alpha.rbegin(), trimmed_alpha.rend());
            // FIXME: depending on the number of orbitals, we should create unsigned long vs dynamic_bitset, so don't use std::stoul per se
            bmqc::SpinString<T> alpha (std::stoul(reversed_alpha, nullptr, 2), this->K);  // convert a binary string to an unsigned long


            // Create a beta SpinString for the second field
            std::string trimmed_beta = boost::algorithm::trim_copy(splitted_line[1]);
            if (trimmed_beta.length() != K) {
                throw std::invalid_argument("ONVExpansion(): One of the provided beta ONVs does not have the correct number of orbitals.");
            }
            std::string reversed_beta(trimmed_beta.rbegin(), trimmed_beta.rend());
            // FIXME: depending on the number of orbitals, we should create unsigned long vs dynamic_bitset, so don't use std::stoul per se
            bmqc::SpinString<T> beta(std::stoul(reversed_beta, nullptr, 2), K);  // convert a binary string to an unsigned long


            // Create a double for the third field
            double coefficient = std::stod(splitted_line[2]);


            // Add an ONVExpansionComponent to the end of this
            this->expansion.emplace_back(ci::ONVExpansionComponent<T>{alpha, beta, coefficient});
        }  // while getline

        this->dim = this->expansion.size();
    }


    /**
     *  Constructor based on a converged @param FCI calculation
     */
    explicit ONVExpansion(const ci::FCI& fci) {

        // Check if the FCI calculation has converged
        if (!fci.is_solved()) {
            throw std::invalid_argument("ONVExpansion(): The given FCI instance isn't solved yet.");
        }

        this->K = fci.get_K();
        this->dim = fci.get_dim();

        Eigen::VectorXd eigenvector = fci.get_eigenvector();

        // Emplace-back this with the alpha-ONVs, beta-ONVs and coefficients
        // The convention that is used is that the alpha addresses are major, i.e. the beta addresses are contiguous
        //      I_alpha I_beta = I_alpha * dim_beta + I_beta
        bmqc::SpinString<T> spin_string_alpha (0, fci.get_addressing_scheme_alpha());
        for (size_t I_alpha = 0; I_alpha < fci.get_dim_alpha(); I_alpha++) {
            if (I_alpha > 0) {
                spin_string_alpha.nextPermutation();
            }

            bmqc::SpinString<T> spin_string_beta (0, fci.get_addressing_scheme_beta());
            for (size_t I_beta = 0; I_beta < fci.get_dim_beta(); I_beta++) {
                if (I_beta > 0) {
                    spin_string_beta.nextPermutation();
                }

                size_t compound_address = I_alpha * fci.get_dim_beta() + I_beta;
                double coefficient = eigenvector(compound_address);

                this->expansion.emplace_back(ci::ONVExpansionComponent<T> {spin_string_alpha, spin_string_beta, coefficient});
            }  // I_beta loop
        }  // I_alpha loop
    }


    /**
     *  Constructor based on a converged @param DOCI calculation
     */
    explicit ONVExpansion(const ci::DOCI& doci) {

        // Check if the DOCI calculation has converged
        if (!doci.is_solved()) {
            throw std::invalid_argument("ONVExpansion(): The given DOCI instance isn't solved yet.");
        }

        this->K = doci.get_K();
        this->dim = doci.get_dim();

        Eigen::VectorXd eigenvector = doci.get_eigenvector();


        // Emplace-back this with the alpha-ONVs, beta-ONVs (which are equal to the alpha-ONVs) and coefficients
        bmqc::SpinString<T> spin_string(0, doci.get_addressing_scheme());
        for (size_t I = 0; I < doci.get_dim(); I++) {
            if (I > 0) {
                spin_string.nextPermutation();
            }

            double coefficient = eigenvector(I);

            this->expansion.emplace_back(ci::ONVExpansionComponent<T>{spin_string, spin_string, coefficient});
        }  // loop over all addresses I
    }


    /*
     *  OPERATORS
     */
    ci::ONVExpansionComponent<T> operator[](size_t index) const {
        return this->expansion[index];
    }


    friend std::ostream& operator<<(std::ostream& os, const ci::ONVExpansion<T>& expansion) {
        for (const auto& component : expansion.expansion) {
            os << component << std::endl;
        }
        return os;
    }


    /*
     *  GETTERS
     */
    std::vector<ci::ONVExpansionComponent<T>> get_expansion() const { return this->expansion; }

    size_t get_dim() const { return this->dim; }

    Eigen::MatrixXd get_one_rdm_aa() const {
        if (!this->are_calculated_one_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this->one_rdm_aa;
    }


    Eigen::MatrixXd get_one_rdm_bb() const {
        if (!this->are_calculated_one_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this-> one_rdm_bb;
    }


    Eigen::MatrixXd get_one_rdm() const {
        if (!this->are_calculated_one_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this-> one_rdm;
    }


    Eigen::Tensor<double, 4> get_two_rdm_aaaa() const {
        if (!this->are_calculated_two_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this->two_rdm_aaaa;
    }


    Eigen::Tensor<double, 4> get_two_rdm_aabb() const {
        if (!this->are_calculated_two_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this->two_rdm_aabb;
    }


    Eigen::Tensor<double, 4> get_two_rdm_bbaa() const {
        if (!this->are_calculated_two_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this->two_rdm_bbaa;
    }


    Eigen::Tensor<double, 4> get_two_rdm_bbbb() const {
        if (!this->are_calculated_two_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this->two_rdm_bbbb;
    }


    Eigen::Tensor<double, 4> get_two_rdm() const {
        if (!this->are_calculated_two_rdms) {
            throw std::logic_error("The requested reduced density matrix is not computed yet.");
        }
        return this->two_rdm;
    }


    /*
     *  PUBLIC METHODS
     */
    /**
     *  @return if this is equal to @param other, within a given @param tolerance
     */
    bool isEqual(const ci::ONVExpansion<T>& other, double tolerance=1.0e-12) const {

        // Check first on equal dimensions
        if (this->dim != other.get_dim()) {
            return false;
        }


        // Check if the SpinStrings are the same
        // In the mean time construct the eigenvector coefficients, so we don't have to loop over the dimension twice
        Eigen::VectorXd this_coefficients = Eigen::VectorXd::Zero(this->dim);
        Eigen::VectorXd other_coefficients = Eigen::VectorXd::Zero(this->dim);
        for (size_t i = 0; i < this->dim; i++) {

            if ((this->expansion[i].alpha != other[i].alpha) || (this->expansion[i].beta != other[i].beta)) {
                return false;
            }

            this_coefficients(i) = this->expansion[i].coefficient;
            other_coefficients(i) = other[i].coefficient;
        }


        // Check if the eigenvectors are equal
        // This is the final check, so we can simplify an if-else construct
        return cpputil::linalg::areEqualEigenvectors(this_coefficients, other_coefficients, tolerance);
    }



    /**
     *  Calculate all the 1-RDMs
     */
    void calculate1RDMs() {

        this->one_rdm_aa = Eigen::MatrixXd::Zero(this->K, this->K);
        this->one_rdm_bb = Eigen::MatrixXd::Zero(this->K, this->K);


        for (size_t I = 0; I < this->dim; I++) {  // loop over all addresses (1)

            bmqc::SpinString<T> alpha_I = this->expansion[I].alpha;
            bmqc::SpinString<T> beta_I = this->expansion[I].beta;
            double c_I = this->expansion[I].coefficient;


            // Calculate the diagonal of the 1-RDMs
            for (size_t p = 0; p < this->K; p++) {

                if (alpha_I.isOccupied(p)) {
                    this->one_rdm_aa(p,p) += std::pow(c_I, 2);
                }

                if (beta_I.isOccupied(p)) {
                    this->one_rdm_bb(p,p) += std::pow(c_I, 2);
                }
            }


            // Calculate the off-diagonal elements, by going over all other ONVs
            for (size_t J = I+1; J < this->dim; J++) {

                bmqc::SpinString<T> alpha_J = this->expansion[J].alpha;
                bmqc::SpinString<T> beta_J = this->expansion[J].beta;
                double c_J = this->expansion[J].coefficient;

                // 1 electron excitation in alpha (i.e. 2 differences), 0 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findOccupiedDifferences(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findOccupiedDifferences(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign, and include it in the RDM contribution
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);
                    this->one_rdm_aa(p,q) += sign * c_I * c_J;
                    this->one_rdm_aa(q,p) += sign * c_I * c_J;
                }


                // 1 electron excitation in beta, 0 in alpha
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = beta_I.findOccupiedDifferences(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = beta_J.findOccupiedDifferences(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign, and include it in the RDM contribution
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);
                    this->one_rdm_bb(p,q) += sign * c_I * c_J;
                    this->one_rdm_bb(q,p) += sign * c_I * c_J;
                }

            }  // loop over addresses J > I
        }  // loop over addresses I


        // The total 1-RDM is the sum of the spin components
        this->one_rdm = this->one_rdm_aa + this->one_rdm_bb;
        this->are_calculated_one_rdms = true;
    }

    /**
     *  Calculate all the 2-RDMS
     */
    void calculate2RDMs() {

        this->two_rdm_aaaa = Eigen::Tensor<double, 4> (this->K, this->K, this->K, this->K);
        this->two_rdm_aaaa.setZero();
        this->two_rdm_aabb = Eigen::Tensor<double, 4> (this->K, this->K, this->K, this->K);
        this->two_rdm_aabb.setZero();
        this->two_rdm_bbaa = Eigen::Tensor<double, 4> (this->K, this->K, this->K, this->K);
        this->two_rdm_bbaa.setZero();
        this->two_rdm_bbbb = Eigen::Tensor<double, 4> (this->K, this->K, this->K, this->K);
        this->two_rdm_bbbb.setZero();


        for (size_t I = 0; I < this->dim; I++) {  // loop over all addresses I

            bmqc::SpinString<T> alpha_I = this->expansion[I].alpha;
            bmqc::SpinString<T> beta_I = this->expansion[I].beta;
            double c_I = this->expansion[I].coefficient;


            for (size_t p = 0; p < this->K; p++) {

                // 'Diagonal' elements of the 2-RDM: aaaa and aabb
                if (alpha_I.isOccupied(p)) {
                    for (size_t q = 0; q < this->K; q++) {
                        if (beta_I.isOccupied(q)) {
                            this->two_rdm_aabb(p,p,q,q) += std::pow(c_I, 2);
                        }

                        if (p != q) {  // can't create/annihilate the same orbital twice
                            if (alpha_I.isOccupied(q)) {
                                this->two_rdm_aaaa(p,p,q,q) += std::pow(c_I, 2);
                                this->two_rdm_aaaa(p,q,q,p) -= std::pow(c_I, 2);
                            }
                        }

                    }  // loop over q
                }

                // 'Diagonal' elements of the 2-RDM: bbbb and bbaa
                if (beta_I.isOccupied(p)) {
                    for (size_t q = 0; q < this->K; q++) {
                        if (alpha_I.isOccupied(q)) {
                            this->two_rdm_bbaa(p,p,q,q) += std::pow(c_I, 2);
                        }

                        if (p != q) {  // can't create/annihilate the same orbital twice
                            if (beta_I.isOccupied(q)) {
                                this->two_rdm_bbbb(p,p,q,q) += std::pow(c_I, 2);
                                this->two_rdm_bbbb(p,q,q,p) -= std::pow(c_I, 2);
                            }
                        }
                    }  // loop over q
                }
            }  // loop over q


            for (size_t J = I+1; J < this->dim; J++) {

                bmqc::SpinString<T> alpha_J = this->expansion[J].alpha;
                bmqc::SpinString<T> beta_J = this->expansion[J].beta;
                double c_J = this->expansion[J].coefficient;


                // 1 electron excitation in alpha, 0 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findOccupiedDifferences(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findOccupiedDifferences(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);


                    for (size_t r = 0; r < this->K; r++) {  // r loops over spatial orbitals

                        if (alpha_I.isOccupied(r) && alpha_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital
                                // Fill in the 2-RDM contributions
                                this->two_rdm_aaaa(p,q,r,r) += sign * c_I * c_J;
                                this->two_rdm_aaaa(r,q,p,r) -= sign * c_I * c_J;
                                this->two_rdm_aaaa(p,r,r,q) -= sign * c_I * c_J;
                                this->two_rdm_aaaa(r,r,p,q) += sign * c_I * c_J;

                                this->two_rdm_aaaa(q,p,r,r) += sign * c_I * c_J;
                                this->two_rdm_aaaa(q,r,r,p) -= sign * c_I * c_J;
                                this->two_rdm_aaaa(r,p,q,r) -= sign * c_I * c_J;
                                this->two_rdm_aaaa(r,r,q,p) += sign * c_I * c_J;
                            }
                        }

                        if (beta_I.isOccupied(r)) {  // beta_I == beta_J from the previous if-branch

                            // Fill in the 2-RDM contributions
                            this->two_rdm_aabb(p,q,r,r) += sign * c_I * c_J;
                            this->two_rdm_aabb(q,p,r,r) += sign * c_I * c_J;

                            this->two_rdm_bbaa(r,r,p,q) += sign * c_I * c_J;
                            this->two_rdm_bbaa(r,r,q,p) += sign * c_I * c_J;
                        }
                    }
                }


                // 0 electron excitations in alpha, 1 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = beta_I.findOccupiedDifferences(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = beta_J.findOccupiedDifferences(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);


                    for (size_t r = 0; r < this->K; r++) {  // r loops over spatial orbitals

                        if (beta_I.isOccupied(r) && beta_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital
                                // Fill in the 2-RDM contributions
                                this->two_rdm_bbbb(p,q,r,r) += sign * c_I * c_J;
                                this->two_rdm_bbbb(r,q,p,r) -= sign * c_I * c_J;
                                this->two_rdm_bbbb(p,r,r,q) -= sign * c_I * c_J;
                                this->two_rdm_bbbb(r,r,p,q) += sign * c_I * c_J;

                                this->two_rdm_bbbb(q,p,r,r) += sign * c_I * c_J;
                                this->two_rdm_bbbb(q,r,r,p) -= sign * c_I * c_J;
                                this->two_rdm_bbbb(r,p,q,r) -= sign * c_I * c_J;
                                this->two_rdm_bbbb(r,r,q,p) += sign * c_I * c_J;
                            }
                        }

                        if (alpha_I.isOccupied(r)) {  // alpha_I == alpha_J from the previous if-branch

                            // Fill in the 2-RDM contributions
                            this->two_rdm_bbaa(p,q,r,r) += sign * c_I * c_J;
                            this->two_rdm_bbaa(q,p,r,r) += sign * c_I * c_J;

                            this->two_rdm_aabb(r,r,p,q) += sign * c_I * c_J;
                            this->two_rdm_aabb(r,r,q,p) += sign * c_I * c_J;
                        }
                    }
                }


                // 1 electron excitation in alpha, 1 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findOccupiedDifferences(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findOccupiedDifferences(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    size_t r = beta_I.findOccupiedDifferences(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t s = beta_J.findOccupiedDifferences(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign, and include it in the 2-RDM contribution
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(s);
                    this->two_rdm_aabb(p,q,r,s) += sign * c_I * c_J;
                    this->two_rdm_aabb(q,p,s,r) += sign * c_I * c_J;

                    this->two_rdm_bbaa(r,s,p,q) += sign * c_I * c_J;
                    this->two_rdm_bbaa(s,r,q,p) += sign * c_I * c_J;
                }


                // 2 electron excitations in alpha, 0 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 4) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = alpha_I.findOccupiedDifferences(alpha_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = alpha_J.findOccupiedDifferences(alpha_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];


                    // Calculate the total sign, and include it in the 2-RDM contribution
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_I.operatorPhaseFactor(r) * alpha_J.operatorPhaseFactor(q) * alpha_J.operatorPhaseFactor(s);
                    this->two_rdm_aaaa(p,q,r,s) += sign * c_I * c_J;
                    this->two_rdm_aaaa(p,s,r,q) -= sign * c_I * c_J;
                    this->two_rdm_aaaa(r,q,p,s) -= sign * c_I * c_J;
                    this->two_rdm_aaaa(r,s,p,q) += sign * c_I * c_J;

                    this->two_rdm_aaaa(q,p,s,r) += sign * c_I * c_J;
                    this->two_rdm_aaaa(s,p,q,r) -= sign * c_I * c_J;
                    this->two_rdm_aaaa(q,r,s,p) -= sign * c_I * c_J;
                    this->two_rdm_aaaa(s,r,q,p) += sign * c_I * c_J;
                }


                // 0 electron excitations in alpha, 2 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 4)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = beta_I.findOccupiedDifferences(beta_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = beta_J.findOccupiedDifferences(beta_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];


                    // Calculate the total sign, and include it in the 2-RDM contribution
                    int sign = beta_I.operatorPhaseFactor(p) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(q) * beta_J.operatorPhaseFactor(s);
                    this->two_rdm_bbbb(p,q,r,s) += sign * c_I * c_J;
                    this->two_rdm_bbbb(p,s,r,q) -= sign * c_I * c_J;
                    this->two_rdm_bbbb(r,q,p,s) -= sign * c_I * c_J;
                    this->two_rdm_bbbb(r,s,p,q) += sign * c_I * c_J;

                    this->two_rdm_bbbb(q,p,s,r) += sign * c_I * c_J;
                    this->two_rdm_bbbb(s,p,q,r) -= sign * c_I * c_J;
                    this->two_rdm_bbbb(q,r,s,p) -= sign * c_I * c_J;
                    this->two_rdm_bbbb(s,r,q,p) += sign * c_I * c_J;
                }

            }  // loop over all addresses J > I

        }  // loop over all addresses I

    this->two_rdm = this->two_rdm_aaaa + this->two_rdm_aabb + this->two_rdm_bbaa + this->two_rdm_bbbb;
    this->are_calculated_two_rdms = true;
    }
};


}  // namespace ci



#endif  // CI_ONVEXPANSION_HPP
