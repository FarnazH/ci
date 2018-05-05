
#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>
#include <string>
#include "ci.hpp"
#include "hf.hpp"
#include "libwint.hpp"
#include "numopt.hpp"

namespace po = boost::program_options;

int main (int argc, char ** argv) {

  std::string filename_xyz;
  std::string basis_set;

  po::variables_map vm;

  try {
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "print help messages")
      ("filename_xyz,f", po::value<std::string>(&filename_xyz)->required(), "filename of the xyz file")
      ("basis_set,b", po::value<std::string>(&basis_set)->required(), "basis set")
      ("charge,c", po::value<int>()->default_value(0), "charge")
      ("rdm,r", po::value<int>()->default_value(0), "reduced density matrix");

    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      std::cout << "DOCI" << std::endl << desc << std::endl;
      std::exit(0);
    }
    po::notify(vm);
  } catch (po::error & e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    return 1;
  } catch(...) {
    std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
  }

  int charge = vm["charge"].as<int>();
  libwint::Molecule mol(filename_xyz, charge);
  libwint::AOBasis ao_basis(mol, basis_set);

  ao_basis.calculateIntegrals();

  hf::rhf::RHF rhf(mol, ao_basis, 1.0e-08);
  rhf.solve();
  
  double internuclear_repulsion_energy = mol.calculateInternuclearRepulsionEnergy();
  double rhf_energy = rhf.get_electronic_energy() + internuclear_repulsion_energy;
  std::cout << "RHF energy: " << std::setprecision(15) << rhf_energy << std::endl;

  Eigen::MatrixXd coefficient_matrix = rhf.get_C_canonical();
  libwint::SOBasis so_basis (ao_basis, coefficient_matrix);

  ci::DOCI doci (so_basis, mol);
  doci.solve(numopt::eigenproblem::SolverType::DAVIDSON);

  double doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;
  std::cout << "DOCI energy: " << std::setprecision(15) << doci_energy << std::endl;

  Eigen::MatrixXd one_rdm;
  Eigen::Tensor<double, 4> two_rdm;
  switch(vm["rdm"].as<int>()) {
  case 0: 
    break;
  case 1: 
    doci.calculate1RDMs();
    one_rdm = doci.get_one_rdm();
    for (int i = 0; i < one_rdm.rows(); i++) {
      for (int j = 0; j < one_rdm.cols(); j++) {
	std::cout << i << "\t" << j << "\t" << one_rdm(i,j) << std::endl;
      }
    }
    break;
  case 2:
    doci.calculate2RDMs();
    two_rdm = doci.get_two_rdm();
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
  }

}
