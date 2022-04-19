/** 
 * @file I1.hpp
 * @brief Implementation of the I sparse integraiton operator
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I1_HPP

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
    * @brief Implementation of the I sparse integration operator
    */ 
   class I1: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          * @param rows Number of rows
          * @param cols Number of columns
          */
         I1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~I1();
         
      protected:

      private:
         /**
          * @brief 1st subdiagonal
          */
         ACoeff_t d_1(const ACoeff_t& n) const; 

         /**
          * @brief 1st superdiagonal
          */
         ACoeff_t d1(const ACoeff_t& n) const; 

         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

}
}
}
}

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I1_HPP
