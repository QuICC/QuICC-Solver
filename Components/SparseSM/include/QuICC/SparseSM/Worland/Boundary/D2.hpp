/** 
 * @file D2.hpp
 * @brief Implementation of the boundary value of second derivative for Worland polynomials
 */

#ifndef QUICC_SPARSESM_WORLAND_BOUNDARY_D2_HPP
#define QUICC_SPARSESM_WORLAND_BOUNDARY_D2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/SparseSM/Worland/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   /**
    * @brief Implementation of the boundary value of second derivative for Worland polynomial
    */ 
   class D2: public ICondition
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         D2(const Scalar_t alpha, const Scalar_t dBeta, const int l);

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

   };

} // Boundary
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_BOUNDARY_D2_HPP
