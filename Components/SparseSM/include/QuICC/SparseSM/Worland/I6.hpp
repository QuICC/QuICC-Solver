/**
 * @file I6.hpp
 * @brief Implementation of the full sphere Worland I6 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I6_HPP
#define QUICC_SPARSESM_WORLAND_I6_HPP

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
#include "QuICC/SparseSM/Worland/I6Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I6 sparse operator
    */
   class I6: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I6(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I6();

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<I6Diags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I6_HPP
