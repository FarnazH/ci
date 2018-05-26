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
#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>
#include <string>

#include "ci.hpp"
#include "hf.hpp"
#include "libwint.hpp"
#include "numopt.hpp"


namespace po = boost::program_options;


/**
 *  Main function for the executable
 */
int main (int argc, char** argv) {

    std::string filename_xyz;
    std::string basis_set;

    po::variables_map variables_map;

    /*
     *  Get the user's input
     */
    try {
        po::options_description desc ("Options");
        desc.add_options()
            ("help,h", "print help messages")
            ("filename_xyz,f", po::value<std::string>(&filename_xyz)->required(), "filename of the xyz file")
            ("basis_set,b", po::value<std::string>(&basis_set)->required(), "basis set name")
            ("charge,c", po::value<int>()->default_value(0), "charge")
            ("rdm,r", po::value<int>()->default_value(0), "reduced density matrices, 1 for 1-RDM, 2 for 2-RDM");

        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::cout << "DOCI" << std::endl << desc << std::endl;
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }

    /*
     *  Do the DOCI calculation
     */

    // Set up the AO basis
    int charge = variables_map["charge"].as<int>();

    libwint::Molecule molecule (filename_xyz, charge);

    libwint::AOBasis ao_basis (molecule, basis_set);
    ao_basis.calculateIntegrals();


    // Do the RHF calculation and prepare an SOBasis
    hf::rhf::RHF rhf(molecule, ao_basis, 1.0e-08);
    rhf.solve();

    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();
    double rhf_energy = rhf.get_electronic_energy() + internuclear_repulsion_energy;

    std::cout << "RHF energy: " << std::setprecision(15) << rhf_energy << std::endl;

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);


    // Do the DOCI calculation
    ci::DOCI doci (so_basis, molecule);

    doci.solve(numopt::eigenproblem::SolverType::DAVIDSON);
    double doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;

    std::cout << "DOCI energy: " << std::setprecision(15) << doci_energy << std::endl;


    // Calculate the DOCI 1- or 2-RDM if requested
    switch(variables_map["rdm"].as<int>()) {
        case 0:
            break;

        case 1:
            doci.calculate1RDMs();
            Eigen::MatrixXd one_rdm = doci.get_one_rdm();
            for (int i = 0; i < one_rdm.rows(); i++) {
                for (int j = 0; j < one_rdm.cols(); j++) {
                    std::cout << i << "\t" << j << "\t" << one_rdm(i,j) << std::endl;
                }
            }
            break;

        case 2:
            doci.calculate2RDMs();
            Eigen::Tensor<double, 4> two_rdm = doci.get_two_rdm();
            for (int i = 0; i < two_rdm.dimensions()[0]; i++) {
                for (int j = 0; j < two_rdm.dimensions()[1]; j++) {
                    for (int k = 0; k < two_rdm.dimensions()[2]; k++) {
                        for (int l = 0; l < two_rdm.dimensions()[3]; l++) {
                            std::cout << i << "\t" << j << "\t" << k << "\t" << l << "\t" << two_rdm(i,j,k,l) << std::endl;
                        }
                    }
                }
            }
            break;

        default:
            break;
    }  // RDM switch

}  // main
