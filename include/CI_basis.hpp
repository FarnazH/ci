#ifndef DOCI_CI_BASIS_HPP
#define DOCI_CI_BASIS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <hf.hpp>
#include "utility.hpp"

namespace doci {

class CI_basis {
//private variables
private:
    Eigen::MatrixXd one_ints;


    private:
        // The one-electron integrals
    Eigen::Tensor<double, 4> two_ints;  // The two-electron integrals
    double internuclear_repulsion;  // The internuclear repulsion energy

    size_t K;  // The number of spatial orbitals
    size_t nelec;  // The number of electrons
//public methods
public:
    /** Default constructor
     */
    CI_basis();

    /** Constructor based on a given RHF instance
     */
    CI_basis(hf::rhf::RHF& rhf);

    /** Constructor based on a given filename
     */
    CI_basis(const std::string& filename);

    /**
     * Getters
     */

    double getOne_int(size_t index1, size_t index2) const;

    double getTwo_int(size_t index1, size_t index2, size_t index3, size_t index4) const;

    double getInternuclear_repulsion() const;

    size_t getK() const;

    size_t getNelec() const;


};

} // namespace doci


#endif // DOCI_CI_BASIS_HPP