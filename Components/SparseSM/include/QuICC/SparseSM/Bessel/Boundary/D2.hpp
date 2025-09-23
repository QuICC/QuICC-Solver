/**
 * @file D2.hpp
 * @brief Implementation of the boundary value of second derivative for Bessel
 * polynomials
 */

#ifndef QUICC_SPARSESM_BESSEL_BOUNDARY_D2_HPP
#define QUICC_SPARSESM_BESSEL_BOUNDARY_D2_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Boundary/ICondition.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

/**
 * @brief Implementation of the boundary value of second derivative for Bessel
 * polynomial
 */
class D2 : public ICondition
{
public:
   /**
    * @brief Constructor for specific alpha,beta pair
    *
    * @param l       Harmonic degree l
    */
   D2(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   ~D2() = default;

   /**
    * @brief Compute list of boundary values
    *
    * @param maxN Highest polynomial
    */
   ACoeff_t compute(const int maxN);

private:
};

} // namespace Boundary
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_BOUNDARY_D2_HPP
