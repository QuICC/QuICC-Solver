/**
 * @file I4.hpp
 * @brief Implementation of the full sphere Worland I4 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I4_HPP
#define QUICC_SPARSESM_WORLAND_I4_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/I4Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I4 sparse operator
    */
   class I4: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          * @param q       Truncation q (only consider rows - q equations)
          */
         I4(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q = 0);

         /**
          * @brief Destructor
          */
         virtual ~I4();

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          *
          * @param list List of triplets (row, col, value)
          */
         virtual void buildTriplets(TripletList_t& list) const override;

         /**
          * @brief Build BLAS banded representation of matrix
          *
          * @param bd   Matrix entries in BLAS banded format
          * @param kL   Number of lower diagonals
          * @param kU   Number of upper diagonals
          */
         virtual void buildBanded(internal::Matrix& bd, unsigned int& kL, unsigned int& kU) const override;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<I4Diags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I4_HPP
