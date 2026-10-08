/**
 * @file OrthogonalValue.cpp
 * @brief Source of the implementation of orthogonal Galerkin stencil for value boundary condition
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "DenseSM/Worland/Stencil/OrthogonalValue.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

#include "Eigen/Dense"
#include "Types/Internal/Typedefs.hpp"
#include <iostream>
namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

OrthogonalValue::Scalar_t OrthogonalValue::basisAlpha(const std::size_t nId)
{
   OrthoId id = static_cast<OrthoId>(nId);

   Scalar_t a;
   switch(id)
   {
      case OrthoId::TorSphEnergy:
         a = Scalar_t(2);
         break;
      case OrthoId::PolSphEnergy:
         a = Scalar_t(1);
         break;
      case OrthoId::ScaSphEnergy:
         a = Scalar_t(2);
         break;
   };


   return a;
}

OrthogonalValue::Scalar_t OrthogonalValue::basisDBeta(const std::size_t nId)
{
   Scalar_t b(0.5);

   return b;
}

OrthogonalValue::OrthogonalValue(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const std::size_t nId, const Scalar_t c) :
    IStencilOperator(rows, cols, alpha, dBeta, l, nId, c)
{
}

void OrthogonalValue::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   OrthoId oid = static_cast<OrthoId>(this->mNId);
   if(this->mL == 0 && (oid == OrthoId::TorSphEnergy || oid == OrthoId::PolSphEnergy))
   {
      mat = Internal::Matrix::Identity(rows,this->rows()-1);
   }
   else
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * this->rows() + this->mL);

      Internal::Array igrid, iweights;
      this->computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::Wnl bW(basisAlpha(this->mNId), basisDBeta(this->mNId));
      Internal::Matrix opBwd(igrid.size(), this->rows()-1);
      Internal::Array diag(opBwd.cols());
      for(int i = 0; i < diag.size(); i++)
      {
         diag(i) = 1./this->norm(i, this->mL);
      }
      bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
         Internal::Array(), ev::Set());
      opBwd = (1.0 - igrid.array().pow(2)).matrix().asDiagonal()*(opBwd * diag.asDiagonal());

      Polynomial::Worland::Wnl W;
      Internal::Matrix opFwd(igrid.size(), this->rows());
      W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
         iweights, ev::Set());

      mat = opFwd.transpose() * opBwd;

      // Prune zeros
      int s = 2;
      for(int i = 0; i < mat.cols()-1; i++)
      {
         mat.block(s+i, i, mat.rows() - (s + i), 1).setZero();
      }
   }
}

OrthogonalValue::Scalar_t OrthogonalValue::norm(const int in, const int il) const
{
   OrthoId id = static_cast<OrthoId>(this->mNId);

   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c;
   switch(id)
   {
      case OrthoId::TorSphEnergy:
         c = l*(l + 1);
         break;
      case OrthoId::PolSphEnergy:
         c = 2*l*(l + 1)*(n + 1)*(2*n + 2*l + 3);
         break;
      case OrthoId::ScaSphEnergy:
         c = 1;
         break;
   };

   return Internal::Math::sqrt(c);
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
