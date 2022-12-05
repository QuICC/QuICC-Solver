/**
 * @file Id.hpp
 * @brief Implementation of the full sphere Worland (restricted) identity sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_ID_HPP
#define QUICC_SPARSESM_WORLAND_ID_HPP

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
#include "QuICC/SparseSM/Worland/IdDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland (restricted) identity sparse operator
    */
   class Id: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Id(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~Id() = default;

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const override;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<IdDiags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_ID_HPP
