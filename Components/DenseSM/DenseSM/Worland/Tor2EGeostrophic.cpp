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

Tor2EGeostrophic::Tor2EGeostrophic(const int nN, const int nL, const int nCpu,
   const Scalar_t alpha, const Scalar_t beta, const bool isGenericBasis,
   const Scalar_t alphaB, const Scalar_t betaB, const bool isTriangular,
   const Scalar_t wAlpha, const Scalar_t wDBeta) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mNcpu(nCpu),
    mUgAlpha(alpha),
    mUgBeta(beta),
    mIsGenericBasis(isGenericBasis),
    mAlphaB(alphaB),
    mBetaB(betaB),
    mIsTriangular(isTriangular),
    mAlpha(wAlpha),
    mDBeta(wDBeta)
{}

Tor2EGeostrophic::Tor2EGeostrophic(const int nN, const int nL, const int nCpu,
   const Scalar_t alpha, const Scalar_t beta, const bool isGenericBasis,
   const Scalar_t alphaB, const Scalar_t betaB, const bool isTriangular) :
    Tor2EGeostrophic(nN, nL, nCpu, alpha, beta, isGenericBasis, alphaB, betaB, isTriangular,
          Polynomial::Worland::worland_default_t::ALPHA,
          Polynomial::Worland::worland_default_t::DBETA)
{}

void Tor2EGeostrophic::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   using namespace Internal::Literals;
   const auto& nN = this->mNn;
   const auto& nL = this->mNl;

   // Compute Z integral quadrature
   const auto nZ = details::GeostrophicTools::cylTruncNz(nN, nL, this->mIsTriangular);
   Internal::Array igridZ;
   Internal::Array iweightZ;
   details::GeostrophicTools::computeQuadratureZ(igridZ, iweightZ, nZ);

   // compute S grid for alpha, beta
   const auto nS = details::GeostrophicTools::cylTruncNs(nN, nL, this->mIsTriangular);
   const auto& alpha = this->mUgAlpha;
   const auto& beta = this->mUgBeta;
   Internal::Array igridS;
   Internal::Array iweightS;
   details::GeostrophicTools::computeGridS(igridS, iweightS, nS, alpha, beta);
   iweightS.resize(0);

   // compute grids for alphaB, betaB
   const auto& alphaB = this->mAlphaB;
   const auto& betaB = this->mBetaB;
   Internal::Array igridX_B, iweightS_B;
   Polynomial::Quadrature::JacobiRule jRuleB(alphaB, betaB);
   jRuleB.computeQuadrature(igridX_B, iweightS_B, nS);

   // grid in s for second set
   Internal::Array igridS_B =
      ((igridX_B.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

   // Compute weights for projecting onto second basis
   int nNug;
   if(this->mIsTriangular)
   {
      nNug = details::GeostrophicTools::cylTruncNug(nL, this->mIsTriangular);
   }
   else
   {
      nNug = details::GeostrophicTools::cylTruncNugC(nN, nL);
   }
   Internal::Matrix iweightS_B_proj;
   iweightS_B_proj.resize(nS, nNug);
   Polynomial::Worland::Wnl wnl(0.5_mp, 0.0_mp);
   wnl.compute<Internal::MHDFloat>(iweightS_B_proj, nNug, 1, igridS_B,
      iweightS_B, Polynomial::Worland::Evaluator::Set());

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
      iweightS_B_proj.col(j).array() *= (1.0_mp - igridX_B.array()) * Cj;
   }

   // allocate integration matrices
   mat.resize(nNug, nN * nL);
   mat.setConstant(0.0);

   // compute integration matrices
   const auto nli = details::GeostrophicTools::nlist(nNug-1, nL);
   const auto& nCpu = this->mNcpu;
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
            details::GeostrophicTools::integrateZ(l, nN - 1, iintgz, igridS,
               igridZ, iweightZ, this->mAlpha, this->mDBeta);
            Internal::Matrix tmp = iweightS_B_proj.transpose() * iintgz;
            mat.block(0, l * nN, nNug, nN) = tmp;
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
