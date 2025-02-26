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
#include "DenseSM/Chebyshev/LinearMap/IProjCrossOperator.hpp"
#include "Types/Internal/Typedefs.hpp"
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
