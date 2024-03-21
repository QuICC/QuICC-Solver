/**
 * @file Tor2EGeostrophic.cpp
 * @brief Source of the implementation of the projection operator from the
 * toroidal scalar to the geostrophic basis
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cstdlib>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Tor2EGeostrophic.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

Tor2EGeostrophic::Tor2EGeostrophic(const int nN, const int nL, const int nS,
   const int nZ, const int nNug, const ArrayI& nli, const int nCpu,
   const Scalar_t alpha, const Scalar_t beta, const bool isGenericBasis,
   const Scalar_t alphaB, const Scalar_t betaB, const bool isTriangular) :
    IEmbeddedOperator(nL * nN, nNug, -0.5, -0.5),
    mNn(nN),
    mNl(nL),
    mNs(nS),
    mNz(nZ),
    mNnug(nNug),
    mNlist(nli),
    mNcpu(nCpu),
    mUgAlpha(alpha),
    mUgBeta(beta),
    mIsGenericBasis(isGenericBasis),
    mAlphaB(alphaB),
    mBetaB(betaB),
    mIsTriangular(isTriangular)
{}

void Tor2EGeostrophic::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   using namespace Internal::Literals;
   const auto& nN = this->mNn;
   const auto& nL = this->mNl;
   const auto& nli = this->mNlist;
   const auto& nS = this->mNs;
   const auto& nZ = this->mNz;
   const auto& nNug = this->mNnug;
   const auto& nCpu = this->mNcpu;
   const auto& alpha = this->mUgAlpha;
   const auto& beta = this->mUgBeta;
   const auto& alphaB = this->mAlphaB;
   const auto& betaB = this->mBetaB;

   Internal::Array igridz;
   Internal::Array iweightz;
   if (QuICCEnv().id() == 0)
   {
      DebuggerMacro_showValue("nz is ", 1, nZ);
      DebuggerMacro_showValue("ns is ", 1, nS);
      DebuggerMacro_showValue("nug is ", 1, nNug);
   }
   details::GeostrophicTools::computeQuadratureZ(igridz, iweightz, nZ);

   // compute Gauss-Jacobi quadrature in x
   Internal::Array igridx, ilambda;
   Polynomial::Quadrature::JacobiRule jRuleA(alpha, beta - 1);
   jRuleA.computeQuadrature(igridx, ilambda, nS);

   // grid in s
   Internal::Array igrids =
      ((igridx.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

   // compute Gauss-Jacobi quadrature in x for second alpha,beta pair
   Internal::Array igridxB, ilambdaB;
   Polynomial::Quadrature::JacobiRule jRuleB(alphaB, betaB);
   jRuleB.computeQuadrature(igridxB, ilambdaB, nS);

   // grid in s for second set
   Internal::Array igridsB =
      ((igridxB.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

   // Compute weights for projecting onto second basis
   Internal::Matrix iEweights_proj;
   iEweights_proj.resize(nS, nNug);
   Polynomial::Worland::Wnl wnl(0.5_mp, 0.0_mp);
   wnl.compute<Internal::MHDFloat>(iEweights_proj, nNug, 1, igridsB,
      ilambdaB, Polynomial::Worland::Evaluator::Set());

   Internal::MHDFloat Cj;
   for (int j = 0; j < nNug; j++)
   {
      if (this->mIsGenericBasis)
      {
         // Jiawen's thesis normalization (via Cnab(j)) 
         // is the same up to the factor 2^(a - a_b - 1) \pi ?
         Cj = Internal::Math::pow(2.0_mp, alpha - alphaB - 1.0_mp)*Internal::Math::PI;
      }
      else
      {
         // Jiawen's thesis normalization (via Cn(j)) 
         // is the same up to the factor sqrt(\pi/(4*sqrt(8))) ?
         Cj = Internal::Math::sqrt(Internal::Math::PI/(4_mp*Internal::Math::sqrt(8_mp)));
      }
      iEweights_proj.col(j).array() *= (1.0_mp - igridxB.array()) * Cj;
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
            // Do nothing
         }
         else
         {
            Internal::Matrix iintgz;
            details::GeostrophicTools::integrateZ(l, nN - 1, iintgz, igrids,
               igridz, iweightz);
            Internal::Matrix tmp = iEweights_proj.transpose() * iintgz;
            mat.block(0, l * nN, nNug, nN) = Internal::cast(tmp);
         }
      }
      pid++;
   }
#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, mat.data(), mat.size(), MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD);
#endif // QUICC_MPI
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
