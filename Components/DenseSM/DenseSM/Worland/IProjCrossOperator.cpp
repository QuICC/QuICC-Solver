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
#include <wigxjpf.h>

// Project includes
//
#include "DenseSM/Worland/IProjCrossOperator.hpp"
#include "Types/Math.hpp"
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

MHDFloat IProjCrossOperator::gaunt(const int lA, const int mA, const int lB,
   const int mB, const int lG, const int mG) const
{
   int lmax = std::max(std::max(lA, lB), lG);

   double val3jA;
   double val3jB;

   wig_table_init(2 * lmax, 3);
   wig_temp_init(2 * lmax);

   /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
    * and 2*m.  To be able to handle half-integer arguments.
    */

   val3jA = wig3jj(2 * lA, 2 * lB, 2 * lG, 2 * mA, 2 * mB, -2 * mG);

   val3jB = wig3jj(2 * lA, 2 * lB, 2 * lG, 0, 0, 0);

   wig_temp_free();
   wig_table_free();

   MHDFloat Kabg = 0;
   if (val3jA != 0 && val3jB != 0)
   {
      MHDFloat ca = static_cast<MHDFloat>(2 * lA + 1);
      MHDFloat cb = static_cast<MHDFloat>(2 * lB + 1);
      MHDFloat cg = static_cast<MHDFloat>(2 * lG + 1);
      Kabg = std::sqrt(ca * cb * cg / (4.0 * Math::PI)) * val3jA * val3jB;
   }

   return Kabg;
}

MHDFloat IProjCrossOperator::elsasser(const int lA, const int mA, const int lB,
   const int mB, const int lG, const int mG) const
{
   int lmax = std::max(std::max(lA, lB + 1), lG);

   double val3jA;
   double val3jB;

   wig_table_init(2 * lmax, 3);
   wig_temp_init(2 * lmax);

   /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
    * and 2*m.  To be able to handle half-integer arguments.
    */

   val3jA = wig3jj(2 * lA, 2 * lB, 2 * lG, 2 * mA, 2 * mB, -2 * mG);

   val3jB = wig3jj(2 * lA, 2 * (lB + 1), 2 * lG, 0, 0, 0);

   wig_temp_free();
   wig_table_free();

   MHDFloat Labg = 0;
   if (val3jA != 0 && val3jB != 0)
   {
      MHDFloat ca = static_cast<MHDFloat>(2 * lA + 1);
      MHDFloat cb = static_cast<MHDFloat>(2 * lB + 1);
      MHDFloat cg = static_cast<MHDFloat>(2 * lG + 1);
      MHDFloat clabg2 = static_cast<MHDFloat>(lA + lB + lG + 2);
      MHDFloat clab_g1 = static_cast<MHDFloat>(lA + lB - lG + 1);
      MHDFloat clbg_a1 = static_cast<MHDFloat>(lB + lG - lA + 1);
      MHDFloat clag_b = static_cast<MHDFloat>(lA + lG - lB);
      Labg = -(std::sqrt(ca * cb * cg / (4.0 * Math::PI)) / 2.0) * val3jA *
             val3jB * std::sqrt(clabg2 * clab_g1 * clbg_a1 * clag_b);
   }

   return Labg;
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
