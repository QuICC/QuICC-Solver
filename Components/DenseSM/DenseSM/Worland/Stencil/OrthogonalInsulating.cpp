/**
 * @file OrthogonalInsulating.cpp
 * @brief Source of the implementation of orthogonal Galerkin stencil for poloidal insulating boundary condition
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "DenseSM/Worland/Stencil/OrthogonalInsulating.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

OrthogonalInsulating::Scalar_t OrthogonalInsulating::basisAlpha(const std::size_t nId)
{
   OrthoId id = static_cast<OrthoId>(nId);

   Scalar_t a;
   switch(id)
   {
      case OrthoId::PolSphEnergy:
         a = Scalar_t(0);
         break;
   };


   return a;
}

OrthogonalInsulating::Scalar_t OrthogonalInsulating::basisDBeta(const std::size_t nId)
{
   Scalar_t b(0.5);

   return b;
}

OrthogonalInsulating::OrthogonalInsulating(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const std::size_t nId, const Scalar_t c) :
    IStencilOperator(rows, cols, alpha, dBeta, l, nId, c)
{
}

void OrthogonalInsulating::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   OrthoId oid = static_cast<OrthoId>(this->mNId);
   if(this->mL == 0 && (oid == OrthoId::PolSphEnergy))
   {
      mat = Internal::Matrix::Identity(rows,this->rows()-1);
   }
   else
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * (this->rows() + 2) + this->mL);

      Internal::Array igrid, iweights;
      this->computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::Wnl bW(basisAlpha(this->mNId), basisDBeta(this->mNId));

      Internal::Matrix opBwd(igrid.size(), this->rows()-1);
      if(oid == OrthoId::PolSphEnergy)
      {
         Internal::Matrix tmp(igrid.size(), this->rows());
         bW.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

         // Build linear combination
         for(int i = 0; i < opBwd.cols(); i++)
         {
            opBwd.col(i) = c1(i, this->mL)*tmp.col(i) + c2(i, this->mL)*tmp.col(i+1);
         }
      }

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

OrthogonalInsulating::Scalar_t OrthogonalInsulating::norm(const OrthoId id, const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c;
   switch(id)
   {
      case OrthoId::PolSphEnergy:
         c = l*(l + 1)*(n + 1)*(n + 2)*(2*n + 2*l + 1)*(2*n + 2*l + 3)*(4*n + 2*l + 5);
         break;
   };

   return Internal::Math::sqrt(c);
}

OrthogonalInsulating::Scalar_t OrthogonalInsulating::c1(const int in, const int il) const
{
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t c = -(n + 2)*(2*n + 2*l + 3)/Internal::Math::sqrt(4*n + 2*l + 3);;

   Scalar_t norm = this->norm(OrthoId::PolSphEnergy, in, il);

   return c/norm;
}

OrthogonalInsulating::Scalar_t OrthogonalInsulating::c2(const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = (n + 1)*(2*n + 2*l + 1)/Internal::Math::sqrt(4*n + 2*l + 7);

   Scalar_t norm = this->norm(OrthoId::PolSphEnergy, in, il);

   return c/norm;
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
