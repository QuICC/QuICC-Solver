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
#include "DenseSM/Worland/Stencil/OrthogonalValue.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

OrthogonalValue::OrthogonalValue(const int rows, const int cols,
   const Scalar_t alpha, const Scalar_t dBeta, const int l) :
    IStencilOperator(rows, cols, alpha, dBeta, l)
{
}

void OrthogonalValue::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (2 * this->rows() + this->mL);

   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

#if 1
   Polynomial::Worland::Wnl bW(this->mcAlpha + 2, Polynomial::Worland::worland_default_t::DBETA);
   Internal::Matrix opBwd(igrid.size(), this->rows()-1);
   bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
      Internal::Array(), ev::Set());
   opBwd = (1.0 - igrid.array().pow(2)).matrix().asDiagonal()*opBwd;
#else
   Polynomial::Worland::Wnl bW(this->mcAlpha + 2, Polynomial::Worland::worland_default_t::DBETA);
   Internal::Matrix tmp(igrid.size(), this->rows());
   bW.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
      Internal::Array(), ev::Set());
   Internal::Matrix opBwd(igrid.size(), this->rows()-1);
   // Set first basis function
   opBwd.col(0) = (1.0 - igrid.array().pow(2)) * igrid.array().pow(this->mL);
   for(int i = 1; i < opBwd.cols(); i++)
   {
      opBwd.col(i) = d1(i+1, this->mL - 1)*tmp.col(i+1) + d2(i+1, this->mL - 1)*tmp.col(i) + d3(i+1, this->mL - 1)*tmp.col(i-1);
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

OrthogonalValue::Scalar_t OrthogonalValue::c1(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t c = (2 * l + 1 + 2 * a + 4 * n) * (2 * l + 5 + 2 * a + 2 * n) * n;

   return c;
}

OrthogonalValue::Scalar_t OrthogonalValue::c2(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t c = -(2*a*a + 2*l*a + 4*n*a + 7*a + 5 + 4 * n *l + 2*l + 4*n*n + 6*n)*(2*l + 3 + 2*a + 4*n);

   return c;
}

OrthogonalValue::Scalar_t OrthogonalValue::c3(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t c = (2*n - 1 + 2*l)*(n + 1 + a) * (2*l + 5 + 2*a + 4*n);

   return c;
}

OrthogonalValue::Scalar_t OrthogonalValue::d1(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t num = n*(1+2*l+2*n)*(2+n+a)*(5+2*l+2*n+2*a);
   Scalar_t den = std::pow(3+2*l+2*n+2*a, 2)*(3+2*l+4*n+2*a)*(7+2*l+4*n+2*a);

   return std::sqrt(num/den);
}

OrthogonalValue::Scalar_t OrthogonalValue::d2(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t num = std::pow(4*n*n+2*l*(1+2*n+a)+(1+a)*(5+2*a)+n*(6+4*a),2);
   Scalar_t den = std::pow((3 + 2*l + 2*n + 2*a)*(1 + 2*l + 4*n + 2*a),2);

   return -std::sqrt(num/den);
}

OrthogonalValue::Scalar_t OrthogonalValue::d3(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   Scalar_t num = (-1+n)*(-1+2*l+2*n)*(1+n+a)*std::pow(5+2*l+4*n+2*a,2);
   Scalar_t den = (3+2*l+2*n+2*a)*(-1+2*l+4*n+2*a)*std::pow(1+2*l+4*n+2*a,2)*(3+2*l+4*n+2*a);

   return std::sqrt(num/den);
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
