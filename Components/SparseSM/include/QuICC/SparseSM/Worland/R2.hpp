/** 
 * @file R2.hpp
 * @brief Implementation of the full sphere Worland R2 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_R2_HPP
#define QUICC_SPARSESM_WORLAND_R2_HPP

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
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/R2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland R2 sparse operator
    */ 
   class R2: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         R2(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~R2();
         
      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<R2Diags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_R2_HPP
