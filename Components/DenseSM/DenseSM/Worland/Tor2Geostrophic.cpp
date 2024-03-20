/**
 * @file Tor2Geostrophic.cpp
 * @brief Source of the implementation of the projection operator from the toroidal scalar to the geostrophic basis
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "Tor2Geostrophic.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/Jacobi/Pnab.hpp"
#include "QuICC/Polynomial/Jacobi/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandChebyshevRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandCylEnergyRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"

#include <iostream>
namespace QuICC {

namespace DenseSM {

namespace Worland {

   Tor2Geostrophic::Tor2Geostrophic(const int nN, const int nL, const int nS, const int nZ, const int nNug, const ArrayI& nli, const int nCpu, const Scalar_t ugAlpha, const Scalar_t ugDBeta, const Scalar_t alphaA, const Scalar_t dBetaA, const bool isTriangular)
      : IGeostrophicOperator(ugAlpha, ugDBeta, nL*nN, nNug, -0.5, -0.5, 0), mNn(nN), mNl(nL), mNs(nS), mNz(nZ), mNnug(nNug), mNlist(nli), mNcpu(nCpu), mAlphaA(alphaA), mDBetaA(dBetaA), mIsTriangular(isTriangular)
   {
   }

   void Tor2Geostrophic::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      using namespace Internal::Literals;
      const auto& nN = this->mNn;
      const auto& nL = this->mNl;
      const auto& nli = this->mNlist;
      const auto& nS = this->mNs;
      const auto& nZ = this->mNz;
      const auto& nNug = this->mNnug;
      const auto& nCpu = this->mNcpu;
      const auto& alphaA = this->mAlphaA;
      const auto& betaA = this->mDBetaA;

   Internal::Array igridz;
   Internal::Array iweightz;
   if (QuICCEnv().id() == 0)
   {
      DebuggerMacro_showValue("nz is ", 1, nZ);
      DebuggerMacro_showValue("ns is ", 1, nS);
      DebuggerMacro_showValue("nug is ", 1, nNug);
   }
   this->computeQuadraturez(igridz, iweightz, nZ);

   // compute Gauss-Jacobi quadrature in x
   Internal::Array igridx, ilambda;
   if (this->isUgBasis())
   {
      Polynomial::Quadrature::JacobiRule jRule(this->mcUgAlpha,
         this->mcUgDBeta - 1.0);
      jRule.computeQuadrature(igridx, ilambda, nS);
   }
   else
   {
      Polynomial::Quadrature::JacobiRule jRuleA(alphaA, betaA);
      jRuleA.computeQuadrature(igridx, ilambda, nS);
   }

   // grid in s
   Internal::Array igrids =
      ((igridx.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

   // compute weights for projecting onto geostrophic basis (integration in s)
   Internal::Matrix ipoly;
   ipoly.resize(nS, nNug);
   Polynomial::Jacobi::Pnab pnab;
   if (this->isUgBasis())
   {
      pnab.compute<Internal::MHDFloat>(ipoly, nNug, this->mcUgAlpha,
         this->mcUgDBeta, igridx, Internal::Array(),
         Polynomial::Jacobi::Evaluator::Set());
   }
   else
   {
      pnab.compute<Internal::MHDFloat>(ipoly, nNug, 0.5_mp,
         1.0_mp, igridx, Internal::Array(),
         Polynomial::Jacobi::Evaluator::Set());
   }

   Internal::Matrix iweights_proj;
   iweights_proj.resize(nS, nNug);
   Internal::MHDFloat Cj;
   for (int j = 0; j < nNug; j++)
   {
      if (this->isUgBasis())
      {
         Internal::MHDFloat a = this->mcUgAlpha;
         Internal::MHDFloat b = this->mcUgDBeta;
         Cj = Internal::Math::sqrt(
            2.0_mp * (2.0_mp * j + a + b + 1.0_mp) *
            Internal::Math::exp(
               Internal::Math::lgamma(j + a + b + 1.0_mp) +
               Internal::Math::lgamma(j + 1.0_mp) -
               Internal::Math::lgamma(j + a + 1.0_mp) -
               Internal::Math::lgamma(j + b + 1.0_mp)));
         iweights_proj.col(j) =
            ilambda.array() *
            ((1.0_mp + igridx.array()) / 2.0_mp).sqrt() *
            Internal::Math::pow(2.0_mp, -a - b - 1.0_mp) * Cj *
            ipoly.col(j).array();
      }
      else
      {
         Cj = Internal::Math::sqrt(
            static_cast<Internal::MHDFloat>((2 * j + 3) * (4 * j + 5)) /
            static_cast<Internal::MHDFloat>(8 * (j + 1)) / Internal::Math::PI);
         iweights_proj.col(j) =
            ilambda.array() * (1.0_mp + igridx.array()).sqrt() *
            Internal::Math::PI / 2.0_mp * Cj * ipoly.col(j).array();
      }
   }

   // allocate integration matrices
   mat.resize(nNug, nN * nL);
   mat.setConstant(0.0);

   // compute integration matrices
   int pid = 0;
   for (int l = 1; l < nL; l += 2)
   {
      if (QuICCEnv().id() == (pid % nCpu))
      {
         if (this->mIsTriangular)
         {
            Internal::Matrix iintgz;
            this->intgzWorland(l, nli(l), iintgz, igrids, igridz, iweightz);
            Internal::Matrix tmp = iweights_proj.transpose() * iintgz;
            mat.block(0, l * nN, nNug, nli(l) + 1) =
               Internal::cast(tmp);
         }
         else
         {
            Internal::Matrix iintgz;
            this->intgzWorland(l, nN - 1, iintgz, igrids, igridz, iweightz);
            Internal::Matrix tmp = iweights_proj.transpose() * iintgz;
            mat.block(0, l * nN, nNug, nN) =
               Internal::cast(tmp);
         }
      }
      pid++;
   }
#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, mat.data(),
      mat.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif // QUICC_MPI
   }

void Tor2Geostrophic::computeQuadraturez(Internal::Array& igridz,
   Internal::Array& iweightz, const int nz) const
{
   using namespace Internal::Literals;

   Polynomial::Quadrature::LegendreRule rule;
   rule.computeQuadrature(igridz, iweightz, nz);

   iweightz = 0.5_mp * iweightz;
}

void Tor2Geostrophic::intgzWorland(int l, int n,
   Internal::Matrix& iintgz, Internal::Array& igrids, Internal::Array& igridz,
   Internal::Array& iweightz) const
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
         r = ((Internal::Math::sqrt(1.0_mp - igrids(k) * igrids(k)) *
                 igridz.array())
                 .square() +
              igrids(k) * igrids(k))
                .sqrt()
                .matrix();
         // cos of theta value
         theta = (Internal::Math::sqrt(1.0_mp - igrids(k) * igrids(k)) *
                  igridz.array() * r.array().inverse())
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

         // compute the integral
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

} // Worland
} // DenseSM
} // QuICC
