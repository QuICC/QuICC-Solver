/**
 * @file OrthogonalValueD1.cpp
 * @brief Source of the implementation of orthogonal Galerkin stencil for value and first derivative boundary condition
 */

// System includes
//
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "DenseSM/Worland/Stencil/OrthogonalValueD1.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

OrthogonalValueD1::Scalar_t OrthogonalValueD1::basisAlpha(const std::size_t nId)
{
   OrthoId id = static_cast<OrthoId>(nId);

   Scalar_t a;
   switch(id)
   {
      case OrthoId::TorSphEnergy:
         a = Scalar_t(4);
         break;
      case OrthoId::PolSphEnergy:
         a = Scalar_t(0);
         break;
      case OrthoId::ScaSphEnergy:
         a = Scalar_t(4);
         break;
      case OrthoId::Lapl2SphEnergy:
         a = Scalar_t(2);
         break;
   };


   return a;
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::basisDBeta(const std::size_t nId)
{
   Scalar_t b(0.5);

   return b;
}

OrthogonalValueD1::OrthogonalValueD1(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const std::size_t nId, const Scalar_t c) :
    IStencilOperator(rows, cols, alpha, dBeta, l, nId, c)
{
}

void OrthogonalValueD1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   OrthoId oid = static_cast<OrthoId>(this->mNId);
   if(this->mL == 0 && (oid == OrthoId::TorSphEnergy || oid == OrthoId::PolSphEnergy || oid == OrthoId::Lapl2SphEnergy))
   {
      mat = Internal::Matrix::Identity(rows,this->rows()-2);
   }
   else
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * (this->rows() + 3) + this->mL);

      Internal::Array igrid, iweights;
      this->computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::Wnl bW(basisAlpha(this->mNId), basisDBeta(this->mNId));

      Internal::Matrix opBwd(igrid.size(), this->rows()-2);
      if(oid == OrthoId::TorSphEnergy || oid == OrthoId::ScaSphEnergy || oid == OrthoId::Lapl2SphEnergy)
      {
         Internal::Array diag(opBwd.cols());
         for(int i = 0; i < diag.size(); i++)
         {
            diag(i) = 1./this->norm(oid, i, this->mL);
         }
         bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());
         opBwd = (1.0 - igrid.array().pow(2)).pow(2).matrix().asDiagonal()*(opBwd * diag.asDiagonal());
      }
      else if(oid == OrthoId::PolSphEnergy)
      {
         Internal::Matrix tmp(igrid.size(), this->rows());
         bW.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

         // Build linear combination
         for(int i = 0; i < opBwd.cols(); i++)
         {
            opBwd.col(i) = c1(i, this->mL)*tmp.col(0) + c2(i, this->mL)*tmp.col(i+1) + c3(i, this->mL)*tmp.col(i+2);
         }
      }

      Polynomial::Worland::Wnl W;
      Internal::Matrix opFwd(igrid.size(), this->rows());
      W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
         iweights, ev::Set());

      mat = opFwd.transpose() * opBwd;

      // Prune zeros
      int s = 3;
      for(int i = 0; i < mat.cols()-1; i++)
      {
         mat.block(s+i, i, mat.rows() - (s + i), 1).setZero();
      }
   }
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::norm(const OrthoId id, const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c;
   switch(id)
   {
      case OrthoId::TorSphEnergy:
         c = l*(l + 1);
         break;
      case OrthoId::PolSphEnergy:
         c = (l*(l + 1)*(n + 1)*(n + 2)*Internal::Math::pow(7+2*l+2*n,2)*Internal::Math::pow(5 + 2*l + 2*n,2))/((5 + 2*l + 2*n)*(7 + 2*l + 2*n)*(9 + 2*l + 4*n));
         break;
      case OrthoId::ScaSphEnergy:
         c = 1;
         break;
      case OrthoId::Lapl2SphEnergy:
         c = l*(l+1)*4*(n + 1)*(n + 2)*(5 + 2*l + 2*n)*(3 + 2*l + 2*n);
         break;
   };

   return Internal::Math::sqrt(c);
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::c1(const int in, const int il) const
{
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = 1./Internal::Math::sqrt(2*l + 3);

   Scalar_t norm = this->norm(OrthoId::PolSphEnergy, in, il);

   return c/norm;
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::c2(const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = (-(2+n)*(7+2*l+2*n))/((9+2*l+4*n)*Internal::Math::sqrt(7 + 2*l + 4*n));

   Scalar_t norm = this->norm(OrthoId::PolSphEnergy, in, il);

   return c/norm;
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::c3(const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = ((1+n)*(5+2*l+2*n))/((9+2*l+4*n)*Internal::Math::sqrt(11 + 2*l + 4*n));

   Scalar_t norm = this->norm(OrthoId::PolSphEnergy, in, il);

   return c/norm;
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
