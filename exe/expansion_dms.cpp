#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/program_options.hpp>

#include <cpputil.hpp>
#include <bmqc.hpp>
#include "ci.hpp"


namespace po = boost::program_options;



/**
 *  Templated function to calculate and output the 1- and 2-RDMs given a wave function expansion in the @param input_file
 */
template <typename T>
void calculateAndOutputDMs(const std::string& input_file) {

    // Create the ONVExpansion and calculate the 1- and 2-RDMs
    ci::ONVExpansion<T> expansion (input_file);
    expansion.calculate1RDMs();
    expansion.calculate2RDMs();


    // Create and open the 1-RDM output files
    std::string output_one_rdm_aa_filename = input_file + "_one_rdm_aa.output";
    std::string output_one_rdm_bb_filename = input_file + "_one_rdm_bb.output";
    std::string output_one_rdm_filename = input_file + "_one_rdm.output";

    std::ofstream output_one_rdm_aa_file;
    output_one_rdm_aa_file.open(output_one_rdm_aa_filename, std::fstream::out);
    cpputil::io::print(expansion.get_one_rdm_aa(), output_one_rdm_aa_file);
    output_one_rdm_aa_file.close();

    std::ofstream output_one_rdm_bb_file;
    output_one_rdm_bb_file.open(output_one_rdm_bb_filename, std::fstream::out);
    cpputil::io::print(expansion.get_one_rdm_bb(), output_one_rdm_bb_file);
    output_one_rdm_bb_file.close();

    std::ofstream output_one_rdm_file;
    output_one_rdm_file.open(output_one_rdm_filename, std::fstream::out);
    cpputil::io::print(expansion.get_one_rdm(), output_one_rdm_file);
    output_one_rdm_file.close();


    // Create and open the 2-RDM output files
    std::string output_two_rdm_aaaa_filename = input_file + "_two_rdm_aaaa.output";
    std::string output_two_rdm_aabb_filename = input_file + "_two_rdm_aabb.output";
    std::string output_two_rdm_bbaa_filename = input_file + "_two_rdm_bbaa.output";
    std::string output_two_rdm_bbbb_filename = input_file + "_two_rdm_bbbb.output";
    std::string output_two_rdm_filename = input_file + "_two_rdm.output";

    std::ofstream output_two_rdm_aaaa_file;
    output_two_rdm_aaaa_file.open(output_two_rdm_aaaa_filename, std::fstream::out);
    cpputil::io::print(expansion.get_two_rdm_aaaa(), output_two_rdm_aaaa_file);
    output_two_rdm_aaaa_file.close();

    std::ofstream output_two_rdm_aabb_file;
    output_two_rdm_aabb_file.open(output_two_rdm_aabb_filename, std::fstream::out);
    cpputil::io::print(expansion.get_two_rdm_aabb(), output_two_rdm_aabb_file);
    output_two_rdm_aabb_file.close();

    std::ofstream output_two_rdm_bbaa_file;
    output_two_rdm_bbaa_file.open(output_two_rdm_bbaa_filename, std::fstream::out);
    cpputil::io::print(expansion.get_two_rdm_bbaa(), output_two_rdm_bbaa_file);
    output_two_rdm_bbaa_file.close();

    std::ofstream output_two_rdm_bbbb_file;
    output_two_rdm_bbbb_file.open(output_two_rdm_bbbb_filename, std::fstream::out);
    cpputil::io::print(expansion.get_two_rdm_bbbb(), output_two_rdm_bbbb_file);
    output_two_rdm_bbbb_file.close();

    std::ofstream output_two_rdm_file;
    output_two_rdm_file.open(output_two_rdm_filename, std::fstream::out);
    cpputil::io::print(expansion.get_two_rdm(), output_two_rdm_file);
    output_two_rdm_file.close();
}




/**
 *  Main function for the executable that calculates 1- and 2-DMs based on a given
 */
int main (int argc, char** argv) {

    // Input processing
    std::string input_file;
    bool more_than_64_orbitals = false;

    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
            ("help,h", "print help messages")
            ("input_file,f", po::value<std::string>(&input_file)->required(), "filename of the input file")
            ("orbitals_flag,o", po::bool_switch(&more_than_64_orbitals), "if the number of orbitals is larger than 64");

        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::cout << "ONV expansion 1- and 2-RDMs" << std::endl << desc << std::endl;
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }


    // Actual calculations
    if (!more_than_64_orbitals) {
        calculateAndOutputDMs<unsigned>(input_file);
    } else {
        calculateAndOutputDMs<boost::dynamic_bitset<>>(input_file);
    }
}
