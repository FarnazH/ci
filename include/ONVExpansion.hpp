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
 *  An ONVExpansion is a decorated std::vector of ci::ONVExpansionComponent<T>s.
 */
template <typename T>
class ONVExpansion : public std::vector<ci::ONVExpansionComponent<T>> {
public:

    // CONSTRUCTORS
    /**
     *  Inherit the bracket initializer from std::vector, see also (https://stackoverflow.com/a/32925226/7930415)
     */
    using std::vector<ci::ONVExpansionComponent<T>>::vector;


    /**
     *  Constructor based on a given @param GAMESS_filename
     */
    explicit ONVExpansion(const std::string& GAMESS_filename) {

        // If the filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
        std::ifstream input_file_stream (GAMESS_filename);

        if (!input_file_stream.good()) {
            throw std::runtime_error("The provided GAMESS file is illegible. Maybe you specified a wrong path?");
        }


        // Do the actual parsing

        // Read in dummy lines up until we actually get to the ONVs and coefficients
        std::string line;
        size_t K = 0;
        while (std::getline(input_file_stream, line)) {
            if (line.find("ALPHA") != std::string::npos && line.find("BETA") != std::string::npos && line.find("COEFFICIENT") != std::string::npos) {  // if find returns an index that's different from the 'not-found' index
                std::getline(input_file_stream, line);

                // After the 'ALPHA', 'BETA', 'COEFFICIENTS' line, we can actually get the number of orbitals by
                // counting the number of dashes (-)
                std::vector<std::string> splitted_line;
                boost::split(splitted_line, line, boost::is_any_of("|"));

                K = splitted_line[0].length();
                break;
            }
        }
        if (K == 0) {
            throw std::invalid_argument("Couldn't read a number of orbitals from a number of dashes.");
        }


        // Read in the ONVs and the coefficients by splitting the line on '|', and then trimming
        // In the GAMESS format, the bit strings are ordered reverse lexicographically, so we will have to reverse
        while (std::getline(input_file_stream, line)) {
            std::vector<std::string> splitted_line;  // create a container for the line to be split in
            boost::split(splitted_line, line, boost::is_any_of("|"));  // split on '|'


            // Create an alpha SpinString for the first field
            std::string trimmed_alpha = boost::algorithm::trim_copy(splitted_line[0]);
            if (trimmed_alpha.length() != K) {
                throw std::invalid_argument("One of the provided alpha ONVs does not have the correct number of orbitals.");
            }
            std::string reversed_alpha (trimmed_alpha.rbegin(), trimmed_alpha.rend());
            bmqc::SpinString<T> alpha (std::stoul(reversed_alpha, nullptr, 2), K);


            // Create a beta SpinString for the second field
            std::string trimmed_beta = boost::algorithm::trim_copy(splitted_line[1]);
            if (trimmed_beta.length() != K) {
                throw std::invalid_argument("One of the provided beta ONVs does not have the correct number of orbitals.");
            }
            std::string reversed_beta (trimmed_beta.rbegin(), trimmed_beta.rend());
            bmqc::SpinString<T> beta (std::stoul(reversed_beta, nullptr, 2), K);


            // Create a double for the third field
            double coefficient = std::stod(splitted_line[2]);


            // Add an ONVExpansionComponent to the end of this
            this->emplace_back(ci::ONVExpansionComponent<T> {alpha, beta, coefficient});
        }  // while
    }


    /**
     *  Constructor based on a converged @param FCI calculation
     */
    explicit ONVExpansion(const ci::FCI& fci) {

        // Check if the FCI calculation has converged
        if (!fci.is_solved()) {
            throw std::invalid_argument("ONVExpansion(): The given FCI instance isn't solved yet.");
        }


        Eigen::VectorXd eigenvector = fci.get_eigenvector();


        // Emplace-back this with the alpha-ONVs, beta-ONVs and coefficients
        // The convention that is used is that the alpha addresses are major, i.e. the beta addresses are contiguous
        //      I_alpha I_beta = I_alpha * dim_beta + I_beta
        bmqc::SpinString<T> spin_string_alpha (1, fci.get_K());
        for (size_t I_alpha = 0; I_alpha < fci.get_dim_alpha(); I_alpha++) {
            if (I_alpha > 0) {
                spin_string_alpha.nextPermutation();
            }

            bmqc::SpinString<T> spin_string_beta (1, fci.get_K());
            for (size_t I_beta = 0; I_beta < fci.get_dim_beta(); I_beta++) {
                if (I_beta > 0) {
                    spin_string_beta.nextPermutation();
                }

                size_t compound_address = I_alpha * fci.get_dim_beta() + I_beta;
                double coefficient = eigenvector(compound_address);

                this->emplace_back(ci::ONVExpansionComponent<T> {spin_string_alpha, spin_string_beta, coefficient});
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


        Eigen::VectorXd eigenvector = doci.get_eigenvector();


        // Emplace-back this with the alpha-ONVs, beta-ONVs (which are equal to the alpha-ONVs) and coefficients
        bmqc::SpinString<T> spin_string (1, doci.get_K());
        for (size_t I = 0; I < doci.get_dim(); I++) {
            if (I > 0) {
                spin_string.nextPermutation();
            }

            double coefficient = eigenvector(I);

            this->emplace_back(ci::ONVExpansionComponent<T> {spin_string, spin_string, coefficient});
        }  // loop over all addresses I
    }


    // PUBLIC METHODS
    /**
     *  @return if this is equal to @param other, within a given @param tolerance
     */
    bool isEqual(const ci::ONVExpansion<T>& other, double tolerance=1.0e-12) const {
        // Check first on equal dimensions
        if (this->size() != other.size()) {
            return false;
        }


        // Check if the SpinStrings are in the same order
        // In the mean time construct the eigenvector coefficients, so we don't have to loop over the dimension twice
        Eigen::VectorXd this_coefficients = Eigen::VectorXd::Zero(this->size());
        Eigen::VectorXd other_coefficients = Eigen::VectorXd::Zero(this->size());
        for (size_t i = 0; i < this->size(); i++) {

            if (((*this)[i].alpha != other[i].alpha) || ((*this)[i].beta != other[i].beta)) {
                return false;
            }

            this_coefficients(i) = (*this)[i].coefficient;
            other_coefficients(i) = other[i].coefficient;
        }


        // Check if the eigenvectors are equal
        // This is the final check, so we can simplify an if-else construct
        return cpputil::linalg::areEqualEigenvectors(this_coefficients, other_coefficients, tolerance);
    }
};


}  // namespace ci




#endif  // CI_ONVEXPANSION_HPP
