/**
 * @file IProjCrossOperator.cpp
 * @brief Source of the implementation of the generic projection of cross
 * product A ^ B
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Worland/IProjCrossOperator.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I3.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

IProjCrossOperator::IProjCrossOperator(const int rows, const int cols, const int q,
   const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, const Scalar_t alpha, const Scalar_t dBeta) :
    IWorlandOperator(rows, cols, alpha, dBeta),
    mIsZero(false),
    mIsImaginary(false),
    mQ(q),
    mLout(lOut),
    mMout(mOut),
    mLa(lA),
    mMa(mA),
    mLb(lB),
    mMb(mB),
    mBand(rows,cols)
{}

bool IProjCrossOperator::isZero() const
{
   return this->mIsZero;
}

bool IProjCrossOperator::isImaginary() const
{
   return this->mIsImaginary;
}

void IProjCrossOperator::applyQI(Internal::Matrix& mat, const int l) const
{
   if(this->mQ > 0)
   {
      // Compute QI
      int nN = mat.rows();
      Internal::SparseMatrix matI;
      std::pair<int,int> band;
      if(this->mQ == 1)
      {
         throw std::logic_error("Quasi-inverse of order 1 is not implemented");
      }
      else if(this->mQ == 2)
      {
         band = std::make_pair(1,3);
         SparseSM::Worland::I2 qi(this->rows(), nN, this->mcAlpha, this->mcDBeta, l);
         matI = qi.mpmat();
      }
      else if(this->mQ == 3)
      {
         band = std::make_pair(1,5);
         SparseSM::Worland::I3 qi(this->rows(), nN, this->mcAlpha, this->mcDBeta, l);
         matI = qi.mpmat();
      }
      else if(this->mQ == 4)
      {
         band = std::make_pair(2,6);
         SparseSM::Worland::I4 qi(this->rows(), nN, this->mcAlpha, this->mcDBeta, l);
         matI = qi.mpmat();
      }
      else
      {
         throw std::logic_error("Quasi-inverse order not implemented");
      }

      // Apply QI
      mat = matI*mat;

      // Get bandwidth of product
      band.first += this->mBand.first;
      band.second += this->mBand.second;

      // Set exact zeros
      for(int i = 0; i < mat.rows(); i++)
      {
         // Set lower part to zero
         for(int j = 0; j < i-band.first; j++)
         {
            mat(i,j) = 0;
         }
         // Set upper part to zero
         for(int j = i + band.second + 1; j < mat.cols(); j++)
         {
            mat(i,j) = 0;
         }
      }
   }
}

void IProjCrossOperator::setBand(const int dL, int s)
{
   if(dL >= 0)
   {
      this->mBand.first = dL; 
      this->mBand.second = 0; 
      s -= 2*dL;
   }
   else
   {
      this->mBand.first = 0; 
      this->mBand.second = -dL; 
   }
   if(s > 0)
   {
      this->mBand.first += s/2;
      this->mBand.second += s/2;
   }
};


} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
