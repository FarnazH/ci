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


public:

    // CONSTRUCTORS
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
        bmqc::SpinString<T> spin_string_alpha(1, this->K);
        for (size_t I_alpha = 0; I_alpha < fci.get_dim_alpha(); I_alpha++) {
            if (I_alpha > 0) {
                spin_string_alpha.nextPermutation();
            }

            bmqc::SpinString<T> spin_string_beta(1, fci.get_K());
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
        bmqc::SpinString<T> spin_string(1, this->K);
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


    /*
     *  GETTERS
     */
    size_t get_dim() const { return this->dim; }


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
};

}  // namespace ci




#endif  // CI_ONVEXPANSION_HPP
