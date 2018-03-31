/**
 *  Create FCI reference data for a given h2 geometry @STO-3G.
 *  For a 2-electron system DOCI is exact, so we can use the current DOCI implementation instead of a FCI one.
 */


#include <libwint.hpp>
#include <numopt.hpp>
#include <hf.hpp>

#include <ci.hpp>



int main() {

    // Do a DOCI calculation based on a RHF calculation
    libwint::Molecule h2 ("../../exe/data/h2.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();

    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-08);
    rhf.solve();
    double rhf_energy = rhf.get_electronic_energy() + internuclear_repulsion_energy;
    std::cout << "RHF energy: " << rhf_energy << std::endl;

    Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
    libwint::SOBasis so_basis (ao_basis, coefficient_matrix);

    ci::DOCI doci (so_basis, h2);
    doci.solve(numopt::eigenproblem::SolverType::DENSE);


    // Calculate the total energy
    double doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
    std::cout << "DOCI energy: " << doci_energy << std::endl;
}
