/** 
 * @file I4D1.hpp
 * @brief Implementation of the I^4 D sparse operator
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4D1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4D1_HPP

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
    * @brief Implementation of the cartesian I^4 D sparse operator, with y = ax + b
    */ 
   class I4D1: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~I4D1();
         
      protected:

      private:
         /**
          * @brief 3rd subdiagonal
          */
         ACoeff_t d_3(const ACoeff_t& n) const; 

         /**
          * @brief 1st subdiagonal
          */
         ACoeff_t d_1(const ACoeff_t& n) const; 

         /**
          * @brief 1st superdiagonal
          */
         ACoeff_t d1(const ACoeff_t& n) const; 

         /**
          * @brief 3rd superdiagonal
          */
         ACoeff_t d3(const ACoeff_t& n) const; 

         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

}
}
}
}

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4D1_HPP
