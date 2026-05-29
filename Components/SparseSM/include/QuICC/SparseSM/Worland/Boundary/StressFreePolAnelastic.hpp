/**
 * @file StressFreePolAnelastic.hpp
 * @brief Implementation of the boundary value of r D2 - r D1(Log(rho))D1 for Worland
 * polynomials
 */

#ifndef QUICC_SPARSESM_WORLAND_BOUNDARY_STRESSFREEPOLANELASTIC_HPP
#define QUICC_SPARSESM_WORLAND_BOUNDARY_STRESSFREEPOLANELASTIC_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D1.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/BasicTypes.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

/**
 * @brief Implementation of the boundary value of second derivative for Worland
 * polynomial
 */
class StressFreePolAnelastic : public ICondition
{
public:
   /**
    * @brief Constructor for specific alpha,beta pair
    *
    * @param alpha   Jacobi alpha
    * @param dBeta   Jacobi beta = l + dBeta
    * @param l       Harmonic degree l
    */
   StressFreePolAnelastic(const Scalar_t alpha, const Scalar_t dBeta, const int l, const MHDFloat Fb);

   /**
    * @brief Destructor
    */
   ~StressFreePolAnelastic() = default;

   /**
    * @brief Compute list of boundary values
    *
    * @param maxN Highest polynomial
    */
   ACoeff_t compute(const int maxN);

private:
   /**
    * @brief Boundary value for k = 0, l = l
    */
   Value mBCk0;

   /**
    * @brief Boundary value for k = 1, l = l+1
    */
   Value mBCk1;

   /**
    * @brief Boundary value for k = 2, l = l+2
    */
   Value mBCk2;

   /**
    * @brief Boundary value First derivative (k=0, l=l)
    */
   D1 mD1;

   /**
    * @brief Radial field boundary value
    */
   MHDFloat mFb;

};

} // namespace Boundary
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_BOUNDARY_STRESSFREEPOLANELASTIC_HPP
