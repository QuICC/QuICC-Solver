/**
 * @file GeostrophicTools.cpp
 * @brief Source of the implementation of the base for a geostrophic projection
 * operator
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "GeostrophicTools.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace details {

Internal::MHDFloat GeostrophicTools::Akj(const int k, const int j)
{
   Internal::MHDFloat dk = static_cast<Internal::MHDFloat>(k);
   Internal::MHDFloat dj = static_cast<Internal::MHDFloat>(j);

   Internal::MHDFloat tmp = 0;
   for (int n = k; n < j + 1; n++)
   {
      Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);
      tmp += Internal::Math::pow(MHD_MP(-1), n) *
             Internal::Math::exp(Internal::Math::lgamma(dj + MHD_MP(1)) -
                                 Internal::Math::lgamma(dj - dn + MHD_MP(1))) /
             (MHD_MP(2) * dn + MHD_MP(1)) *
             Internal::Math::exp(
                Internal::Math::lgamma(MHD_MP(2) * dn + MHD_MP(2)) +
                Internal::Math::lgamma(dn + dk + MHD_MP(2)) -
                Internal::Math::lgamma(dn + MHD_MP(1)) -
                Internal::Math::lgamma(dn - dk + MHD_MP(1)) -
                Internal::Math::lgamma(
                   MHD_MP(2) * dk + MHD_MP(2) * dn + MHD_MP(4)));
   }

   return 2 * Internal::Math::PI *
          Internal::Math::sqrt(
             (MHD_MP(4) * dk + MHD_MP(3)) / (MHD_MP(4) * Internal::Math::PI)) *
          Internal::Math::pow(MHD_MP(2), (2 * k + 2)) * tmp;
}

Internal::MHDFloat GeostrophicTools::Bjnab(const int j, const int n,
   const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::MHDFloat dj = static_cast<Internal::MHDFloat>(j);
   Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);

   auto Cn = [&](const Internal::MHDFloat dn)
   {
      // normalisation factor
      return Internal::Math::sqrt(
         (MHD_MP(2) * (MHD_MP(2) * dn + a + b + MHD_MP(1))) *
         Internal::Math::exp(Internal::Math::lgamma(dn + a + b + MHD_MP(1)) +
                             Internal::Math::lgamma(dn + MHD_MP(1)) -
                             Internal::Math::lgamma(dn + a + MHD_MP(1)) -
                             Internal::Math::lgamma(dn + b + MHD_MP(1))));
   };

   Internal::MHDFloat tmp = 0;
   for (int m = j; m < n + 1; m++)
   {
      Internal::MHDFloat dm = static_cast<Internal::MHDFloat>(m);
      Internal::MHDFloat bFactor =
         Internal::Math::exp(Internal::Math::lgamma(dn + MHD_MP(2)) +
                             Internal::Math::lgamma(dm + MHD_MP(2)) -
                             Internal::Math::lgamma(dn - dm + MHD_MP(1)) -
                             Internal::Math::lgamma(dm + MHD_MP(1)) -
                             Internal::Math::lgamma(dm - dj + MHD_MP(1)) -
                             Internal::Math::lgamma(dj + MHD_MP(1)));
      tmp += MHD_MP(1) / ((dn + 1) * (dm + 1)) * bFactor *
             Internal::Math::exp(
                Internal::Math::lgamma(a + b + dn + dm + MHD_MP(1)) -
                Internal::Math::lgamma(a + dm + MHD_MP(1))) *
             Internal::Math::pow(-MHD_MP(1), m - j);
   }

   Internal::MHDFloat ret =
      Cn(dn) *
      Internal::Math::exp(Internal::Math::lgamma(a + dn + MHD_MP(1)) -
                          Internal::Math::lgamma(dn + MHD_MP(1)) -
                          Internal::Math::lgamma(a + b + dn + MHD_MP(1))) *
      tmp;
   return ret;
}

Internal::MHDFloat GeostrophicTools::Bjn(const int j, const int n)
{
   Internal::MHDFloat a = MHD_MP(0.5);
   Internal::MHDFloat b = MHD_MP(1);
   Internal::MHDFloat dj = static_cast<Internal::MHDFloat>(j);
   Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);

   Internal::MHDFloat tmp = 0;
   for (int m = j; m < n + 1; m++)
   {
      Internal::MHDFloat dm = static_cast<Internal::MHDFloat>(m);
      Internal::MHDFloat bFactor =
         Internal::Math::exp(Internal::Math::lgamma(dn + MHD_MP(2)) +
                             Internal::Math::lgamma(dm + MHD_MP(2)) -
                             Internal::Math::lgamma(dn - dm + MHD_MP(1)) -
                             Internal::Math::lgamma(dm + MHD_MP(1)) -
                             Internal::Math::lgamma(dm - dj + MHD_MP(1)) -
                             Internal::Math::lgamma(dj + MHD_MP(1)));
      tmp += MHD_MP(1) / ((dn + 1) * (dm + 1)) * bFactor *
             Internal::Math::exp(
                Internal::Math::lgamma(a + b + dn + dm + MHD_MP(1)) -
                Internal::Math::lgamma(a + dm + MHD_MP(1))) *
             Internal::Math::pow(-MHD_MP(1), m - j);
   }

   Internal::MHDFloat ret =
      Internal::Math::sqrt(
         (MHD_MP(2) * dn + MHD_MP(3)) * (MHD_MP(4) * dn + MHD_MP(5)) /
         (MHD_MP(8) * Internal::Math::PI * (dn + MHD_MP(1)))) *
      Internal::Math::exp(Internal::Math::lgamma(a + dn + MHD_MP(1)) -
                          Internal::Math::lgamma(dn + MHD_MP(1)) -
                          Internal::Math::lgamma(a + b + dn + MHD_MP(1))) *
      tmp;
   return ret;
}

Internal::MHDFloat GeostrophicTools::Cnab(const int n,
   const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   using namespace Internal::Literals;

   Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);
   Internal::MHDFloat ret =
      Internal::Math::sqrt((2_mp * (2_mp * dn + a + b + 1_mp))) *
      Internal::Math::exp(0.5_mp * (Internal::Math::lgamma(dn + a + b + 1_mp) +
                                      Internal::Math::lgamma(dn + 1_mp) -
                                      Internal::Math::lgamma(dn + a + 1_mp) -
                                      Internal::Math::lgamma(dn + b + 1_mp)));
   return ret;
}

Internal::MHDFloat GeostrophicTools::Cn(const int n)
{
   using namespace Internal::Literals;

   Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);
   Internal::MHDFloat ret =
      Internal::Math::sqrt((2_mp * dn + 3_mp) * (4_mp * dn + 5_mp) /
                           (8_mp * Internal::Math::PI * (dn + 1_mp)));
   return ret;
}

void GeostrophicTools::computeQuadratureZ(Internal::Array& igridz,
   Internal::Array& iweightz, const int nz)
{
   using namespace Internal::Literals;

   Polynomial::Quadrature::LegendreRule rule;
   rule.computeQuadrature(igridz, iweightz, nz);

   iweightz = 0.5_mp * iweightz;
}

void GeostrophicTools::integrateZ(int l, int n, Internal::Matrix& iintgz,
   Internal::Array& igrids, Internal::Array& igridz, Internal::Array& iweightz)
{
   using namespace Internal::Literals;

   assert(l / 2 + n + 1 <= igridz.size());
   int nz = igridz.size();
   int ns = igrids.size();
   iintgz.resize(ns, n + 1);

   if (l % 2 == 1)
   {
      Internal::Array r;
      Internal::Array theta;
      r.resize(nz);
      theta.resize(nz);

      // Loop over geostrophic cylinders
      for (int k = 0; k < igrids.size(); k++)
      {
         const auto& s_ = igrids(k);
         // z in [-1, 1] mapped to sphere at s: \hat{z}(s) = sqrt(1-s^2) z
         // r = sqrt(\hat{z}^2 + s^2) = sqrt((sqrt(1 - s^2)*z)^2 + s^2)
         r = ((Internal::Math::sqrt(1.0_mp - s_ * s_) * igridz.array())
                 .square() +
              s_ * s_)
                .sqrt()
                .matrix();
         // cos of theta value: cos\theta = \hat{z}/r = sqrt(1-s^2)z/r
         theta = (Internal::Math::sqrt(1.0_mp - s_ * s_) * igridz.array() *
                  r.array().inverse())
                    .matrix();

         // compute worland values on cylinder
         Internal::Matrix ipoly;
         ipoly.resize(nz, n + 1);
         Polynomial::Worland::Wnl wnl;
         wnl.compute<Internal::MHDFloat>(ipoly, n + 1, l, r, Internal::Array(),
            Polynomial::Worland::Evaluator::Set());

         // compute derivative of legendre poly on cylinder
         Internal::Matrix ipolyP;
         Internal::Matrix idiff;
         ipolyP.resize(nz, l + 1);
         idiff.resize(nz, l + 1);
         Polynomial::ALegendre::Plm plm;
         plm.compute<Internal::MHDFloat>(ipolyP, l + 1, 0, theta,
            Internal::Array(), Polynomial::ALegendre::Evaluator::Set());
         Polynomial::ALegendre::dPlm dPlm;
         dPlm.compute<Internal::MHDFloat>(idiff, l + 1, 0, theta,
            Internal::Array(), Polynomial::ALegendre::Evaluator::Set());

         // compute the integral: v_\phi = -d_\theta T_l^0
         Internal::Array tmp =
            -(iweightz.array() * idiff.col(l).array()).matrix();
         iintgz.row(k) = tmp.transpose() * ipoly;
      }
   }
   else
   {
      iintgz.setConstant(0.0_mp);
   }
}

int GeostrophicTools::cylTruncNug(const int nL, const bool isTriangular)
{
   int nN;
   if (isTriangular)
   {
      nN = int(((nL - 2) - (nL - 2) % 2) / 2);
   }
   else
   {
      nN = int(((nL - 2) - (nL - 2) % 2) / 2);
   }

   return nN;
}

int GeostrophicTools::cylTruncNugC(const int nr, const int nL)
{
   int nN = int(((nL - 2) - (nL - 2) % 2) / 2) + nr;

   return nN;
}

int GeostrophicTools::cylTruncNs(const int nr, const int nL,
   const bool isTriangular)
{
   int nN;
   if (isTriangular)
   {
      nN = int(((nL + 2) - (nL + 2) % 2) / 2 + 3) + 1;
   }
   else
   {
      nN = int(((nL - 1) - (nL - 1) % 2) / 2 + nr + 1) + 1;
   }

   return nN;
}

int GeostrophicTools::cylTruncNz(const int nr, const int nL,
   const bool isTriangular)
{
   int nN;
   if (isTriangular)
   {
      nN = int(((nL - 1) - (nL - 1) % 2) / 2 + 2) + 1;
   }
   else
   {
      nN = int(((nL - 1) - (nL - 1) % 2) / 2 + nr + 1) + 1;
   }

   return nN;
}

int GeostrophicTools::cylTruncNr(const int nL, const bool isTriangular)
{
   int nN;
   if (isTriangular)
   {
      nN = int(((nL + 1) - (nL + 1) % 2) / 2 + 2) + 1;
   }
   else
   {
      nN = int(nL + 1 + 2);
   }

   return nN;
}

ArrayI GeostrophicTools::nlist(const int maxNug, const int maxnl)
{
   if (2 * maxNug + 1 > maxnl - 1)
   {
      throw std::logic_error("L truncation is not enough to capture all "
                             "geostrophic modes required");
   }

   ArrayI nli = ArrayI::Constant(maxnl, 0);

   for (int l = 0; l < maxnl; l++)
   {
      if (l % 2 == 0)
      {
         nli[l] = -1;
      }
   }

   for (int k = 0; k <= maxNug; k++)
   {
      nli(2 * k + 1) = maxNug - k;
   }

   return nli;
}

} // namespace details
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
