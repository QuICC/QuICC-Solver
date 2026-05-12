/**
 * @file OrthogonalValueD1.cpp
 * @brief Source of the implementation of orthogonal Galerkin stencil for value and first derivative boundary condition
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/Stencil/OrthogonalValueD1.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

OrthogonalValueD1::OrthogonalValueD1(const int rows, const int cols,
   const Scalar_t alpha, const Scalar_t dBeta, const int l) :
    IStencilOperator(rows, cols, alpha, dBeta, l)
{
}

void OrthogonalValueD1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (2 * (this->rows() + 3) + this->mL);

   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

#if 1
   Polynomial::Worland::Wnl bW(this->mcAlpha + 4, Polynomial::Worland::worland_default_t::DBETA);
   Internal::Matrix opBwd(igrid.size(), this->rows()-2);
   bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
      Internal::Array(), ev::Set());
   opBwd = (1.0 - igrid.array().pow(2)).pow(2).matrix().asDiagonal()*opBwd;
#else
   Polynomial::Worland::Wnl bW(this->mcAlpha + 2, Polynomial::Worland::worland_default_t::DBETA);
   Internal::Matrix tmp(igrid.size(), this->rows() + 3);
   bW.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
      Internal::Array(), ev::Set());
   Internal::Matrix opBwd(igrid.size(), this->rows()-2);
   // Build linear combination
   for(int i = 0; i < opBwd.cols(); i++)
   {
      opBwd.col(i) = d1(i+1, this->mL - 1)*tmp.col(i+2) + d2(i+1, this->mL - 1)*tmp.col(i+1) + d3(i+1, this->mL - 1)*tmp.col(i);
   }
#endif

   Polynomial::Worland::Wnl W;
   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
      iweights, ev::Set());

   mat = opFwd.transpose() * opBwd;
   Internal::Array norm = ((mat.transpose() * mat)).diagonal().array().sqrt().pow(-1).matrix();

   mat = mat * norm.asDiagonal();
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::c1(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t c = (n + 1)*n*(2*l + 5 + 2*a + 4*n);

   return c;
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::c2(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t c = -2*n*(n + a + 3)*(2*a + 4*n + 7 + 2*l);

   return c;
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::c3(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t c = (n + a + 3) * (n + 2 + a) * (2*l + 9 + 2*a + 4*n);

   return c;
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::d1(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t num = n*(1+n)*(3+2*l+2*n)*std::pow(5+2*l+4*n+2*a,2);
   Scalar_t den = 4*(1+2*l+2*n)*std::pow(7+2*l+2*n+2*a,2)*(7+2*l+4*n+2*a)*(11+2*l+4*n+2*a);

   return std::sqrt(num/den);
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::d2(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t num = n*(3 + n + a);
   Scalar_t den = (1 + 2*l + 2*n)*(7 + 2*l + 2*n + 2*a);

   return -std::sqrt(num/den);
}

OrthogonalValueD1::Scalar_t OrthogonalValueD1::d3(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t num = (2+n+a)*(3+n+a)*(5+2*l+2*n+2*a)*std::pow(9+2*l+4*n+2*a,2);
   Scalar_t den = 4*std::pow(1+2*l+2*n,2)*(7+2*l+2*n+2*a)*(3+2*l+4*n+2*a)*(7+2*l+4*n+2*a);

   return std::sqrt(num/den);
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
