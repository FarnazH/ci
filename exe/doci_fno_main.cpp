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


    // Prepare an initial SOBasis by using Löwdin orthogonalization
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes_S (ao_basis.get_S());
    libwint::SOBasis so_basis (ao_basis, saes_S.operatorInverseSqrt());  // Löwdin orthogonalization of the AOBasis

    // Get the FCI natural orbitals
    ci::FCI fci (so_basis, molecule.get_N()/2, molecule.get_N()/2);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

    fci.calculate1RDMs();
    Eigen::MatrixXd D = fci.get_one_rdm();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes_D (D);

    Eigen::MatrixXd U = saes_D.eigenvectors();
    so_basis.rotate(U);


    // Do a DOCI calculation based on the FCI natural orbital basis
    ci::DOCI doci (so_basis, molecule);
    doci.solve(&dense_options);


    // Print the energy to an output file
    // Create and open a file: filename.xyz -> filename_doci_rhf_basisset.output
    std::string output_filename = xyz_filename;
    boost::replace_last(output_filename, ".xyz", std::string("_doci_fno_") + basisset + std::string(".output"));

    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);

    output_file << std::setprecision(15) << doci.get_eigenvalue() + internuclear_repulsion_energy << std::endl;

    output_file.close();
}
