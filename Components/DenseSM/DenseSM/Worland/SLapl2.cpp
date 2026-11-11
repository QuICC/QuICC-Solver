/**
 * @file SLapl2.cpp
 * @brief Source of the implementation of spherical bilaplacian tau operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/SLapl2.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

SLapl2::Scalar_t SLapl2::basisAlpha(const std::size_t nId)
{
   OrthoId id = static_cast<OrthoId>(nId);

   Scalar_t a;
   switch(id)
   {
      case OrthoId::PolNoSlipSphEnergy:
      case OrthoId::PolNoSlipSphEnergyWorland:
         a = Scalar_t(1);
         break;
      case OrthoId::PolInsulatingSphEnergy:
      case OrthoId::PolInsulatingSphEnergyWorland:
         a = Scalar_t(0);
         break;
      case OrthoId::Worland:
         a = Scalar_t(-2.0);
         break;
   };


   return a;
}

SLapl2::Scalar_t SLapl2::basisDBeta(const std::size_t nId)
{
   OrthoId id = static_cast<OrthoId>(nId);

   Scalar_t b;
   if(id == OrthoId::Worland)
   {
      b = Scalar_t(-2.0);
   }
   else
   {
      b = Scalar_t(0.5);
   }

   return b;
}

SLapl2::SLapl2(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const std::size_t bId) :
    ITauOperator(rows, cols, alpha, dBeta, l, bId)
{
}

SLapl2::SLapl2(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const OrthoId bId) :
    ITauOperator(rows, cols, alpha, dBeta, l, static_cast<std::size_t>(bId))
{
}

void SLapl2::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   OrthoId bid = static_cast<OrthoId>(this->mBId);

   if(bid == OrthoId::Worland)
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * this->rows() + this->mL);

      Internal::Array igrid, iweights;
      this->computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::slaplWnl bW;
      Internal::Matrix opBwd(igrid.size(), this->rows());
      bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

      Polynomial::Worland::Wnl W;
      Internal::Matrix opFwd(igrid.size(), this->rows());
      W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
            iweights, ev::Set());

      mat = opFwd.transpose() * opBwd;

      // Prune zeros
      for(int i = 0; i < mat.cols(); i++)
      {
         mat.block(i, i, mat.rows() - i, 1).setZero();
      }

      mat = mat*mat;
   }
   else
   {
      throw std::logic_error("Not yet implemented");
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
