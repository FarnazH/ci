#ifndef CI_ONVEXPANSIONCOMPONENT_HPP
#define CI_ONVEXPANSIONCOMPONENT_HPP


#include <iostream>

#include <bmqc.hpp>


namespace ci {

/**
 *  A templated struct that holds the linear expansion @param coefficient, together with the SpinString<T>
 *  representation of the ONVs, separated in @param alpha and @param beta components
 */
template <typename T>
struct ONVExpansionComponent {
    bmqc::SpinString<T> alpha;
    bmqc::SpinString<T> beta;
    double coefficient;


    bool isEqual(const ci::ONVExpansionComponent<T>& other, double tolerance = 1.0e-12) {
        return (this->alpha == other.alpha)
               && (this->beta == other.beta)
               && (std::abs(this->coefficient - other.coefficient) < tolerance);
    }


    friend std::ostream& operator<<(std::ostream& os, const ci::ONVExpansionComponent<T>& expansion_component) {
        os << expansion_component.alpha << '|' << expansion_component.beta << ' ' << expansion_component.coefficient;
        return os;
    }
};


}  // namespace ci


#endif  // CI_ONVEXPANSIONCOMPONENT_HPP
