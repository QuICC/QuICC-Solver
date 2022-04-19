/** 
 * @file I4D1R1.hpp
 * @brief Implementation of the full sphere Worland I4D1R1 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I4D1R1_HPP
#define QUICC_SPARSESM_WORLAND_I4D1R1_HPP

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
#include "QuICC/SparseSM/Worland/I4D1R1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I4D1R1 sparse operator
    */ 
   class I4D1R1: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4D1R1(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I4D1R1();
         
      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<I4D1R1Diags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I4D1R1_HPP
