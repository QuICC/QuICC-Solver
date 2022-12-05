/** 
 * @file Value.hpp
 * @brief Implementation of the bounary value for Worland polynomials
 */

#ifndef QUICC_SPARSESM_WORLAND_BOUNDARY_VALUE_HPP
#define QUICC_SPARSESM_WORLAND_BOUNDARY_VALUE_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/SparseSM//Worland/I2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   /**
    * @brief Implementation of the Worland polynomial
    */ 
   class Value: public IDiags
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          */
         Value(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         ~Value() = default;

         /**
          * @brief Compute list of boundary values
          *
          * @param maxN Highest polynomial 
          */
         ACoeff_t compute(const int maxN, const int k = 0, const bool normalized = true);

      private:

   };

} // Boundary
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_BOUNDARY_VALUE_HPP
