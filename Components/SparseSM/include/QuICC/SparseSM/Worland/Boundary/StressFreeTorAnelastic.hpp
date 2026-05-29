/**
 * @file StressFreeTorAnelastic.hpp
 * @brief Implementation of the boundary value of r D 1/r - D1(Log(rho)) for Worland
 * polynomials
 */

#ifndef QUICC_SPARSESM_WORLAND_BOUNDARY_STRESSFREETORANELASTIC_HPP
#define QUICC_SPARSESM_WORLAND_BOUNDARY_STRESSFREETORANELASTIC_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/BasicTypes.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

/**
 * @brief Implementation of the boundary value of r D 1/r - D1(Log(rho)) for Worland polynomial
 */
class StressFreeTorAnelastic : public ICondition
{
public:
   /**
    * @brief Constructor for specific alpha,beta pair
    *
    * @param alpha   Jacobi alpha
    * @param dBeta   Jacobi beta = l + dBeta
    * @param l       Harmonic degree l
    */
   StressFreeTorAnelastic(const Scalar_t alpha, const Scalar_t dBeta, const int l, const MHDFloat Fb);

   /**
    * @brief Destructor
    */
   ~StressFreeTorAnelastic() = default;

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
    * @brief Radial field boundary value
    */
   MHDFloat mFb;
};

} // namespace Boundary
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_BOUNDARY_STRESSFREETORANELASTIC_HPP
