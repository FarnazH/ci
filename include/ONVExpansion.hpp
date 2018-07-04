#ifndef CI_ONVEXPANSION_HPP
#define CI_ONVEXPANSION_HPP


#include "ONVExpansionComponent.hpp"

#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/range/adaptors.hpp>


namespace ci {


/**
 *  An ONVExpansion is a decorated std::vector of ci::ONVExpansionComponent<T>s.
 */
template <typename T>
class ONVExpansion : public std::vector<ci::ONVExpansionComponent<T>> {
public:

    // CONSTRUCTORS
    /**
     *  Inherit the bracket initializer, see also (https://stackoverflow.com/a/32925226/7930415)
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
};


}  // namespace ci




#endif  // CI_ONVEXPANSION_HPP
