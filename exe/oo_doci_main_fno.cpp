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


    // Get the FCI natural orbitals
    ci::FCI fci (so_basis, molecule.get_N()/2, molecule.get_N()/2);
    // Specify solver options and solve the eigenvalue problem
    numopt::eigenproblem::DenseSolverOptions dense_options;
    fci.solve(&dense_options);

    fci.calculate1RDMs();
    Eigen::MatrixXd D = fci.get_one_rdm();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (D);

    Eigen::MatrixXd U = saes.eigenvectors();
    so_basis.rotate(U);


    // Do a DOCI calculation
    ci::DOCI doci (so_basis, molecule);
    doci.solve(&dense_options);


    // Specify solver options and perform the orbital optimization
    numopt::eigenproblem::DavidsonSolverOptions davidson_options;
    //  In lexical notation, the Hartree-Fock determinant has the lowest address
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(doci.get_dim());

    initial_guess(0) = 1;
    davidson_options.X_0 = initial_guess;
    doci.orbitalOptimize(&davidson_options);


    // Print the energy to an output file
    // Create and open a file: filename.xyz -> filename_fci.output
    std::string output_filename = xyz_filename;
    boost::replace_last(output_filename, ".xyz", std::string("_oo_doci_fno_") + basisset + std::string(".output"));

    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);

    output_file << std::setprecision(15) << doci.get_eigenvalue() + internuclear_repulsion_energy << std::endl;

    output_file.close();
}
