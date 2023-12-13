/**
 * @file IEmbeddedSMOperator.cpp
 * @brief Source of the implementation of generic interface to the dense
 * operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/IEmbeddedSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

IEmbeddedSMOperator::IEmbeddedSMOperator(const int rows, const int cols) :
    mRows(rows), mCols(cols)
{}

int IEmbeddedSMOperator::rows() const
{
   return this->mRows;
}

int IEmbeddedSMOperator::cols() const
{
   return this->mCols;
}

Matrix IEmbeddedSMOperator::mat() const
{
   Internal::Matrix mat;
   this->buildOpImpl(mat, this->rows(), this->cols());

   return mat.cast<MHDFloat>();
}

Matrix IEmbeddedSMOperator::embedded(const int r, const int c) const
{
   Internal::Matrix mat;
   this->buildOpImpl(mat, r, c);

   return mat.cast<MHDFloat>();
}

} // namespace DenseSM
} // namespace QuICC
