/**
 * @file I4Lapl.hpp
 * @brief Implementation of the full sphere Worland I4Lapl sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I4LAPL_HPP
#define QUICC_SPARSESM_WORLAND_I4LAPL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/I4LaplDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I4Lapl sparse operator
    */
   class I4Lapl: public IWorlandOperator
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
         I4Lapl(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q = 0);

         /**
          * @brief Destructor
          */
         virtual ~I4Lapl() = default;

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          *
          * @param list List of triplets (row, col, value)
          */
         void buildTriplets(TripletList_t& list) const final;

         /**
          * @brief Build BLAS banded representation of matrix
          *
          * @param bd   Matrix entries in BLAS banded format
          * @param kL   Number of lower diagonals
          * @param kU   Number of upper diagonals
          */
         void buildBanded(internal::Matrix& bd, unsigned int& kL, unsigned int& kU) const final;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<I4LaplDiags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I4LAPL_HPP
