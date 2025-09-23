/**
 * @file Value.hpp
 * @brief Implementation of the boundary value for Bessel polynomials
 */

#ifndef QUICC_SPARSESM_BESSEL_BOUNDARY_VALUE_HPP
#define QUICC_SPARSESM_BESSEL_BOUNDARY_VALUE_HPP

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
 * @brief Implementation of the boundary value for Bessel polynomial
 */
class Value : public ICondition
{
public:
   /**
    * @brief Constructor for specific alpha,beta pair
    *
    * @param l       Harmonic degree l
    */
   Value(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   ~Value() = default;

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

#endif // QUICC_SPARSESM_BESSEL_BOUNDARY_VALUE_HPP
