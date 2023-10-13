/**
 * @file IDenseSMOperator.cpp
 * @brief Source of the implementation of generic interface to the dense operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/IDenseSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

   IDenseSMOperator::IDenseSMOperator(const int rows, const int cols)
      : mRows(rows), mCols(cols)
   {
   }

   int IDenseSMOperator::rows() const
   {
      return this->mRows;
   }

   int IDenseSMOperator::cols() const
   {
      return this->mCols;
   }

   Matrix IDenseSMOperator::mat() const
   {
      Internal::Matrix mat;
      this->buildOpImpl(mat, this->rows(), this->cols());

      return mat.cast<MHDFloat>();
   }

   Matrix IDenseSMOperator::embedded(const int r, const int c) const
   {
      Internal::Matrix mat;
      this->buildOpImpl(mat, r, c);

      return mat.cast<MHDFloat>();
   }

}
}
