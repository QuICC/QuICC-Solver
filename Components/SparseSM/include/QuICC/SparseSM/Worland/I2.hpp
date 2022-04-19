/**
 * @file I2.hpp
 * @brief Implementation of the full sphere Worland I2 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I2_HPP
#define QUICC_SPARSESM_WORLAND_I2_HPP

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
#include "QuICC/SparseSM/Worland/I2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I2 sparse operator
    */
   class I2: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I2();

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const override;

         /**
          * @brief Build BLAS banded representation of matrix
          */
         virtual void buildBanded(internal::Matrix& bd, unsigned int& kL, unsigned int& kU) const override;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<I2Diags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I2_HPP
