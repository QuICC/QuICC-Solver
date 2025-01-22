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
#include "DenseSM/Chebyshev/LinearMap/IProjCrossOperator.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Math.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

IProjCrossOperator::IProjCrossOperator(const int rows, const int cols,
   const int q, const int p, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, const Scalar_t lower, const Scalar_t upper) :
    ILinearMapOperator(rows, cols, lower, upper),
    mIsZero(false),
    mIsImaginary(false),
    mQ(q), mP(p),
    mLout(lOut),
    mMout(mOut),
    mLa(lA),
    mMa(mA),
    mLb(lB),
    mMb(mB)
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

void IProjCrossOperator::applyQI(Internal::Matrix& mat) const
{
   if(this->mQ > 0)
   {
      int nN = mat.rows();
      Internal::SparseMatrix matI;
      if(this->mQ == 1)
      {
         SparseSM::Chebyshev::LinearMap::I1 qi(this->rows(), nN, this->mcLower, this->mcUpper);
         matI = qi.mpmat();
      }
      else if(this->mQ == 2)
      {
         SparseSM::Chebyshev::LinearMap::I2 qi(this->rows(), nN, this->mcLower, this->mcUpper);
         matI = qi.mpmat();
      }
      else if(this->mQ == 3)
      {
         SparseSM::Chebyshev::LinearMap::I3 qi(this->rows(), nN, this->mcLower, this->mcUpper);
         matI = qi.mpmat();
      }
      else if(this->mQ == 4)
      {
         SparseSM::Chebyshev::LinearMap::I4 qi(this->rows(), nN, this->mcLower, this->mcUpper);
         matI = qi.mpmat();
      }
      else
      {
         throw std::logic_error("Quasi-inverse order not implemented");
      }

      mat = matI*mat;
   }
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
