#include <fstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>

#include "ci.hpp"
#include "hf.hpp"
#include "libwint.hpp"
#include "numopt.hpp"


int main (int argc, char** argv) {

    // Input processing
    // Argument 1: .xyz-filename
    std::string xyz_filename = argv[1];

    // Argument 2: basisset
    std::string basisset = argv[2];


    // Actual calculations
    // Prepare the AO basis
    libwint::Molecule molecule (xyz_filename);
    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();

    libwint::AOBasis ao_basis (molecule, basisset);
    ao_basis.calculateIntegrals();


    // Prepare the SO basis from RHF coefficients
    hf::rhf::RHF rhf (molecule, ao_basis, 1.0e-06);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    libwint::SOBasis so_basis (ao_basis, rhf.get_C_canonical());


    // Do a FCI calculation
    ci::FCI fci (so_basis, molecule.get_N()/2, molecule.get_N()/2);  // assume N_alpha = N_beta
    fci.solve(numopt::eigenproblem::SolverType::DENSE);


    // Print the energy to an output file
    // Create and open a file: filename.xyz -> filename_fci_lowdin_basisset.output
    std::string output_filename = xyz_filename;
    boost::replace_last(output_filename, ".xyz", std::string("_fci_rhf_") + basisset + std::string(".output"));

    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);

    output_file << std::setprecision(15) << fci.get_eigenvalue() + internuclear_repulsion_energy << std::endl;

    output_file.close();
}
