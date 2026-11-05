/**
 * @file Id.cpp
 * @brief Source of the implementation of dense identity operator operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "DenseSM/Worland/Id.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

Id::Scalar_t Id::basisAlpha(const std::size_t nId)
{
   OrthoId id = static_cast<OrthoId>(nId);

   Scalar_t a;
   switch(id)
   {
      case OrthoId::TorValueSphEnergy:
      case OrthoId::TorValueSphEnergyWorland:
         a = Scalar_t(2);
         break;
      case OrthoId::PolNoSlipSphEnergy:
      case OrthoId::PolNoSlipSphEnergyWorland:
      case OrthoId::PolNoSlipSphEnergyNoLaplWorland:
         a = Scalar_t(0);
         break;
      case OrthoId::PolInsulatingSphEnergy:
      case OrthoId::PolInsulatingSphEnergyWorland:
         a = Scalar_t(0);
         break;
      case OrthoId::ScaValueSphEnergy:
      case OrthoId::ScaValueSphEnergyWorland:
         a = Scalar_t(2);
         break;
      case OrthoId::Worland:
         a = Scalar_t(-2.0);
         break;
   };


   return a;
}

Id::Scalar_t Id::basisDBeta(const std::size_t nId)
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

Id::Id(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const std::size_t bId) :
    ITauOperator(rows, cols, alpha, dBeta, l, bId)
{
}

void Id::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   OrthoId bid = static_cast<OrthoId>(this->mBId);

   if(bid == OrthoId::Worland ||
      bid == OrthoId::ScaValueSphEnergy ||
      bid == OrthoId::TorValueSphEnergy||
      bid == OrthoId::PolInsulatingSphEnergy ||
      bid == OrthoId::PolNoSlipSphEnergy)
   {
      mat = Internal::Matrix::Identity(rows, cols);
   }
   else if(bid == OrthoId::ScaValueSphEnergyWorland)
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * (this->rows()+3) + this->mL);

      Internal::Array igrid, iweights;
      Polynomial::Worland::worland_sphenergy_t::Rule quad;
      quad.computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::Wnl bW;
      Internal::Matrix opBwd(igrid.size(), this->rows());
      bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

      Polynomial::Worland::Wnl W(basisAlpha(this->mBId), basisDBeta(this->mBId));
      Internal::Matrix opFwd(igrid.size(), this->rows());
      Internal::Array diag(opBwd.cols());
      for(int i = 0; i < diag.size(); i++)
      {
         diag(i) = 1./this->norm(bid, i, this->mL);
      }
      W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
         iweights, ev::Set());
      opFwd = (1.0 - igrid.array().pow(2)).matrix().asDiagonal()*(opFwd * diag.asDiagonal());

      mat = opFwd.transpose() * opBwd;
   }
   else if(bid == OrthoId::TorValueSphEnergyWorland)
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * (this->rows()+3) + this->mL);

      Internal::Array igrid, iweights;
      Polynomial::Worland::worland_sphenergy_t::Rule quad;
      quad.computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::Wnl bW;
      Internal::Matrix opBwd(igrid.size(), this->rows());
      bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

      Polynomial::Worland::Wnl W(basisAlpha(this->mBId), basisDBeta(this->mBId));
      Internal::Matrix opFwd(igrid.size(), this->rows());
      Internal::Array diag(opBwd.cols());
      for(int i = 0; i < diag.size(); i++)
      {
         diag(i) = 1./this->norm(bid, i, this->mL);
      }
      W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
         iweights, ev::Set());
      opFwd = (1.0 - igrid.array().pow(2)).matrix().asDiagonal()*(opFwd * diag.asDiagonal());

      Scalar_t laplh = Scalar_t(this->mL*(this->mL+1));
      mat = laplh*(opFwd.transpose() * opBwd);
   }
   else if(bid == OrthoId::PolNoSlipSphEnergyWorland)
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * (this->rows()+3) + this->mL);

      Internal::Array igrid, iweights;
      Polynomial::Worland::worland_sphenergy_t::Rule quad;
      quad.computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::slaplWnl bW;
      Internal::Matrix opBwd(igrid.size(), this->rows());
      bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

      Polynomial::Worland::Wnl W(basisAlpha(this->mBId), basisDBeta(this->mBId));
      Internal::Matrix opFwd(igrid.size(), this->rows());
      Internal::Matrix tmp(igrid.size(), this->rows()+2);
      W.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
         iweights, ev::Set());

      // Build linear combination
      for(int i = 0; i < opFwd.cols(); i++)
      {
         opFwd.col(i) = c1(i, this->mL)*tmp.col(0) + c2(i, this->mL)*tmp.col(i+1) + c3(i, this->mL)*tmp.col(i+2);
      }

      Scalar_t laplh = -Scalar_t(this->mL*(this->mL+1));
      mat = laplh*(opFwd.transpose() * opBwd);
   }
   else if(bid == OrthoId::PolNoSlipSphEnergyNoLaplWorland)
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * (this->rows()+3) + this->mL);

      Internal::Array igrid, iweights;
      Polynomial::Worland::worland_sphenergy_t::Rule quad;
      quad.computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::Wnl bW;
      Internal::Matrix opBwd(igrid.size(), this->rows());
      bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

      Polynomial::Worland::Wnl W(basisAlpha(this->mBId), basisDBeta(this->mBId));
      Internal::Matrix opFwd(igrid.size(), this->rows());
      Internal::Matrix tmp(igrid.size(), this->rows()+2);
      W.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
         iweights, ev::Set());

      // Build linear combination
      for(int i = 0; i < opFwd.cols(); i++)
      {
         opFwd.col(i) = c1(i, this->mL)*tmp.col(0) + c2(i, this->mL)*tmp.col(i+1) + c3(i, this->mL)*tmp.col(i+2);
      }

      Scalar_t laplh = -Scalar_t(this->mL*(this->mL+1));
      mat = laplh*(opFwd.transpose() * opBwd);
   }
   else if(bid == OrthoId::PolInsulatingSphEnergyWorland)
   {
      namespace ev = Polynomial::Worland::Evaluator;
      const int nR = (2 * (this->rows()+2) + this->mL);

      Internal::Array igrid, iweights;
      Polynomial::Worland::worland_sphenergy_t::Rule quad;
      quad.computeQuadrature(igrid, iweights, nR);

      Polynomial::Worland::slaplWnl bW;
      Internal::Matrix opBwd(igrid.size(), this->rows());
      bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
            Internal::Array(), ev::Set());

      Internal::Matrix opFwd(igrid.size(), this->rows());
      Internal::Matrix tmp(igrid.size(), this->rows()+1);
      Polynomial::Worland::Wnl W(basisAlpha(this->mBId), basisDBeta(this->mBId));
      W.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
         iweights, ev::Set());

      // Build linear combination
      for(int i = 0; i < opFwd.cols(); i++)
      {
         opFwd.col(i) = d1(i, this->mL)*tmp.col(i) + d2(i, this->mL)*tmp.col(i+1);
      }

      Scalar_t laplh = Scalar_t(this->mL*(this->mL+1));
      mat = -laplh*(opFwd.transpose() * opBwd);

      // Add boundary term
      Internal::Array bgrid = Internal::Array::Ones(1);
      tmp.resize(bgrid.size(), this->rows()+1);
      W.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, bgrid,
         Internal::Array(), ev::Set());
      for(int i = 0; i < opFwd.cols(); i++)
      {
         opFwd.col(i) = d1(i, this->mL)*tmp.col(i) + d2(i, this->mL)*tmp.col(i+1);
      }

      opBwd.resize(bgrid.size(), this->rows());
      Polynomial::Worland::Wnl bcW;
      bcW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, bgrid,
            Internal::Array(), ev::Set());
      Polynomial::Worland::dWnl bcdW;
      for(int j = 0; j < mat.cols(); j++)
      {
         for(int i = 0; i < mat.rows(); i++)
         {
            mat(i,j) += laplh*Scalar_t(this->mL+1)*opFwd(0,i)*opBwd(0,j);
         }
      }
      bcdW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, bgrid,
            Internal::Array(), ev::Set());
      for(int j = 0; j < mat.cols(); j++)
      {
         for(int i = 0; i < mat.rows(); i++)
         {
            mat(i,j) += laplh*opFwd(0,i)*opBwd(0,j);
         }
      }
   }
   else
   {
      throw std::logic_error("Not yet implemented");
   }
}

Id::Scalar_t Id::norm(const OrthoId id, const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c;
   switch(id)
   {
      case OrthoId::Worland:
         c = 1.0;
         break;
      case OrthoId::TorValueSphEnergy:
      case OrthoId::TorValueSphEnergyWorland:
         c = l*(l + 1);
         break;
      case OrthoId::ScaValueSphEnergy:
      case OrthoId::ScaValueSphEnergyWorland:
         c = 1.0;
         break;
      case OrthoId::PolInsulatingSphEnergy:
      case OrthoId::PolInsulatingSphEnergyWorland:
         c = l*(l + 1)*(n + 1)*(n + 2)*(2*n + 2*l + 1)*(2*n + 2*l + 3)*(4*n + 2*l + 5);
         break;
      case OrthoId::PolNoSlipSphEnergy:
      case OrthoId::PolNoSlipSphEnergyWorland:
      case OrthoId::PolNoSlipSphEnergyNoLaplWorland:
         c = (l*(l + 1)*(n + 1)*(n + 2)*Internal::Math::pow(7+2*l+2*n,2)*Internal::Math::pow(5 + 2*l + 2*n,2))/((5 + 2*l + 2*n)*(7 + 2*l + 2*n)*(9 + 2*l + 4*n));
         break;
   };

   return Internal::Math::sqrt(c);
}

Id::Scalar_t Id::c1(const int in, const int il) const
{
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = 1./Internal::Math::sqrt(2*l + 3);

   Scalar_t norm = this->norm(OrthoId::PolNoSlipSphEnergy, in, il);

   return c/norm;
}

Id::Scalar_t Id::c2(const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = (-(2+n)*(7+2*l+2*n))/((9+2*l+4*n)*Internal::Math::sqrt(7 + 2*l + 4*n));

   Scalar_t norm = this->norm(OrthoId::PolNoSlipSphEnergy, in, il);

   return c/norm;
}

Id::Scalar_t Id::c3(const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = ((1+n)*(5+2*l+2*n))/((9+2*l+4*n)*Internal::Math::sqrt(11 + 2*l + 4*n));

   Scalar_t norm = this->norm(OrthoId::PolNoSlipSphEnergy, in, il);

   return c/norm;
}

Id::Scalar_t Id::d1(const int in, const int il) const
{
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t c = -(n + 2)*(2*n + 2*l + 3)/Internal::Math::sqrt(4*n + 2*l + 3);;

   Scalar_t norm = this->norm(OrthoId::PolInsulatingSphEnergy, in, il);

   return c/norm;
}

Id::Scalar_t Id::d2(const int in, const int il) const
{
   Scalar_t n = static_cast<Scalar_t>(in);
   Scalar_t l = static_cast<Scalar_t>(il);
   Scalar_t c = (n + 1)*(2*n + 2*l + 1)/Internal::Math::sqrt(4*n + 2*l + 7);

   Scalar_t norm = this->norm(OrthoId::PolInsulatingSphEnergy, in, il);

   return c/norm;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
