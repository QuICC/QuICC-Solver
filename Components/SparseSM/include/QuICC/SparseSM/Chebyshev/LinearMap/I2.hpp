/** 
 * @file I2.hpp
 * @brief Implementation of the I^2 sparse integration operator
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the I^2 sparse integration operator
    */ 
   class I2: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~I2();
         
      protected:

      private:
         /**
          * @brief 2nd subdiagonal
          */
         ACoeff_t d_2(const ACoeff_t& n) const; 

         /**
          * @brief diagonal
          */
         ACoeff_t d0(const ACoeff_t& n) const; 

         /**
          * @brief 2nd superdiagonal
          */
         ACoeff_t d2(const ACoeff_t& n) const; 

         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

}
}
}
}

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2_HPP
