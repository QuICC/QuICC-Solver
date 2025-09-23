/**
 * @file D1.hpp
 * @brief Implementation of the boundary value of first derivative for Bessel
 * polynomials
 */

#ifndef QUICC_SPARSESM_BESSEL_BOUNDARY_D1_HPP
#define QUICC_SPARSESM_BESSEL_BOUNDARY_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Boundary/ICondition.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

/// Namespace for Boundary tau lines for Bessel polynomials
namespace Boundary {

/**
 * @brief Implementation of the boundary value of first derivative for Bessel
 * polynomial
 */
class D1 : public ICondition
{
public:
   /**
    * @brief Constructor for specific alpha,beta pair
    *
    * @param l       Harmonic degree l
    */
   D1(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   ~D1() = default;

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

#endif // QUICC_SPARSESM_BESSEL_BOUNDARY_D1_HPP
