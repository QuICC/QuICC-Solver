/**
 * @file IMatrixSMOperator.cpp
 * @brief Source of the implementation of generic interface to the dense
 * spectral operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/IMatrixSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

IMatrixSMOperator::IMatrixSMOperator(const int rows, const int cols) :
    mRows(rows), mCols(cols)
{}

int IMatrixSMOperator::rows() const
{
   return this->mRows;
}

int IMatrixSMOperator::cols() const
{
   return this->mCols;
}

Matrix IMatrixSMOperator::mat() const
{
   Internal::Matrix mat;
   this->buildOpImpl(mat, this->rows(), this->cols());

   return mat.cast<MHDFloat>();
}

} // namespace DenseSM
} // namespace QuICC
