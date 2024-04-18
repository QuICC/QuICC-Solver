/**
 * @file Tor2Geostrophic.cpp
 * @brief Source of the implementation of the projection operator from the
 * toroidal scalar to the geostrophic basis
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "Tor2Geostrophic.hpp"
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

Tor2Geostrophic::Tor2Geostrophic(const int nN, const int nL, const int nCpu,
   const Scalar_t alpha, const Scalar_t beta, const bool isGenericBasis,
   const bool isTriangular, const Scalar_t wAlpha, const Scalar_t wDBeta) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mNcpu(nCpu),
    mUgAlpha(alpha),
    mUgBeta(beta),
    mIsGenericBasis(isGenericBasis),
    mIsTriangular(isTriangular),
    mAlpha(wAlpha),
    mDBeta(wDBeta)
{}

Tor2Geostrophic::Tor2Geostrophic(const int nN, const int nL, const int nCpu,
   const Scalar_t alpha, const Scalar_t beta, const bool isGenericBasis,
   const bool isTriangular) :
    Tor2Geostrophic(nN, nL, nCpu, alpha, beta, isGenericBasis, isTriangular, 
          Polynomial::Worland::worland_default_t::ALPHA,
          Polynomial::Worland::worland_default_t::DBETA)
{}

void Tor2Geostrophic::buildOpImpl(Internal::Matrix& mat, const int rows,
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

   // Compute S grid
   const auto nS = details::GeostrophicTools::cylTruncNs(nN, nL, this->mIsTriangular);
   const auto& alpha = this->mUgAlpha;
   const auto& beta = this->mUgBeta;
   Internal::Array igridS;
   Internal::Array iweightS;
   details::GeostrophicTools::computeGridS(igridS, iweightS, nS, alpha, beta);

   // Compute integrator for geostrophic basis (weighted projector)
   int nNug;
   if(this->mIsTriangular)
   {
      nNug = details::GeostrophicTools::cylTruncNug(nL, this->mIsTriangular);
   }
   else
   {
      nNug = details::GeostrophicTools::cylTruncNugC(nN, nL);
   }
   Internal::Matrix iweightS_proj;
   iweightS_proj.resize(nS, nNug);
   Internal::MHDFloat c;
   if (this->mIsGenericBasis)
   {
      // Jiawen's thesis normalization (via Cnab(j)) 
      // is the same up to the factor 2^-(a + 2) ?
      c = Internal::Math::pow(2.0_mp, -alpha - 2.0_mp);
   }
   else
   {
      // Jiawen's thesis normalization (via Cn(j)) 
      // is the same up to the factor sqrt(\pi/8) ?
      c = Internal::Math::sqrt(Internal::Math::PI/8_mp);
   }
   Polynomial::Worland::Wnl wnl(alpha, beta-1);
   wnl.compute<Internal::MHDFloat>(iweightS_proj, nNug, 1, igridS,
      iweightS * c, Polynomial::Worland::Evaluator::Set());

   // allocate integration matrices
   mat.resize(nNug, nN * nL);
   mat.setConstant(0.0);

   // compute integration matrices
   const auto nli = details::GeostrophicTools::nlist(nNug - 1, nL);
   const auto& nCpu = this->mNcpu;
   int pid = 0;
   for (int l = 1; l < nL; l += 2)
   {
      if (QuICCEnv().id() == (pid % nCpu))
      {
         int maxN_;
         if (this->mIsTriangular)
         {
            maxN_ = nli(l);
         }
         else
         {
            maxN_ = nN - 1;
         }

         Internal::Matrix iintgz;
         details::GeostrophicTools::integrateZ(l, maxN_, iintgz, igridS, igridZ,
            iweightZ, this->mAlpha, this->mDBeta);
         Internal::Matrix tmp = iweightS_proj.transpose() * iintgz;
         mat.block(0, l * nN, nNug, maxN_ + 1) = tmp;
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
