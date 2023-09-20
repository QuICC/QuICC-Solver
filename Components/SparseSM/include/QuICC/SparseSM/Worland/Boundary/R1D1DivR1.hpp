/**
 * @file R1D1DivR1.hpp
 * @brief Implementation of the boundary value of r D 1/r for Worland polynomials
 */

#ifndef QUICC_SPARSESM_WORLAND_BOUNDARY_R1D1DIVR1_HPP
#define QUICC_SPARSESM_WORLAND_BOUNDARY_R1D1DIVR1_HPP

// System includes
//

// Project includes
//
#include "Types/Precision.hpp"
#include "QuICC/SparseSM/Worland/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   /**
    * @brief Implementation of the boundary value of r D 1/r for Worland polynomial
    */
   class R1D1DivR1: public ICondition
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         R1D1DivR1(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         ~R1D1DivR1() = default;

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

   };

} // Boundary
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_BOUNDARY_R1D1DIVR1_HPP
