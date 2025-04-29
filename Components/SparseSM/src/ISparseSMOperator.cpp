/**
 * @file ISparseSMOperator.cpp
 * @brief Source of the implementation of generic interface to the sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/ISparseSMOperator.hpp"

namespace QuICC {

namespace SparseSM {

   ISparseSMOperator::ISparseSMOperator(const int rows, const int cols)
      : mRows(rows), mCols(cols)
   {
   }

   int ISparseSMOperator::rows() const
   {
      return this->mRows;
   }

   int ISparseSMOperator::cols() const
   {
      return this->mCols;
   }

   SparseMatrix ISparseSMOperator::mat() const
   {
      Internal::SparseMatrix mat;
      this->buildOpImpl(mat, this->rows(), this->cols());

      return mat.cast<MHDFloat>();
   }

   SparseMatrix ISparseSMOperator::embedded(const int r, const int c) const
   {
      Internal::SparseMatrix mat;
      this->buildOpImpl(mat, r, c);

      return mat.cast<MHDFloat>();
   }

   Matrix ISparseSMOperator::banded(unsigned int& kL, unsigned int &kU) const
   {
      Internal::Matrix bd;
      this->buildBanded(bd, kL, kU);

      return bd.cast<MHDFloat>();
   }

   void ISparseSMOperator::buildOpImpl(Internal::SparseMatrix& mat, const int r, const int c) const
   {
      // Build triplets
      TripletList_t list;
      this->buildTriplets(list);

      // Build sparse matrix
      mat.resize(r, c);
      mat.setFromTriplets(list.begin(), list.end());
      mat.makeCompressed();
   }

   void ISparseSMOperator::buildOpImpl(Internal::Matrix& bd, const int r, const int c) const
   {
      unsigned int kL, kU;
      this->buildBanded(bd, kL, kU);
      if(kU > 0)
      {
         bd.topLeftCorner(1,1)(0,0) = kU;
      }
      if(kL > 0)
      {
         bd.bottomRightCorner(1,1)(0,0) = kL;
      }
   }

   void ISparseSMOperator::buildBanded(Internal::Matrix& bd, unsigned int& kL, unsigned int &kU) const
   {
      throw std::logic_error("Banded operator not yet implemented!");
   }

   void ISparseSMOperator::convertToTriplets(TripletList_t& list, const int d, const ACoeffI& rowIdx, const ACoeff_t& diag) const
   {
      for(int i = 0; i < diag.size(); ++i)
      {
         int row = rowIdx(i);
         int col = rowIdx(i) + d;
         if(col < 0)
         {
            this->leftOutOfMatrix(list, row, col, diag(i));
         } else if(col >= this->cols())
         {
            this->rightOutOfMatrix(list, row, col, diag(i));
         } else
         {
            if(diag(i) != MHD_MP(0.0))
            {
               list.push_back(Triplet_t(row, col, diag(i)));
            }
         }
      }
   }

   void ISparseSMOperator::leftOutOfMatrix(TripletList_t& , const int , const int , const Scalar_t ) const
   {
      // Ignore entry
   }

   void ISparseSMOperator::rightOutOfMatrix(TripletList_t& , const int , const int , const Scalar_t ) const
   {
      // Ignore entry
   }

}
}
