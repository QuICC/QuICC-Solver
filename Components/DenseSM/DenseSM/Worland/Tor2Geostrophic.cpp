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
#include "QuICC/QuICCEnv.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

Tor2Geostrophic::Tor2Geostrophic(const int nN, const int nL, const int nCpu,
   const Scalar_t alpha, const Scalar_t beta, const bool isGenericBasis,
   const bool isTriangular) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mNcpu(nCpu),
    mUgAlpha(alpha),
    mUgBeta(beta),
    mIsGenericBasis(isGenericBasis),
    mIsTriangular(isTriangular)
{}

void Tor2Geostrophic::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   using namespace Internal::Literals;
   const auto& nN = this->mNn;
   const auto& nL = this->mNl;
   const auto nS = details::GeostrophicTools::cylTruncNs(nN, nL, this->mIsTriangular);
   const auto nZ = details::GeostrophicTools::cylTruncNz(nN, nL, this->mIsTriangular);
   int nNug;
   if(this->mIsTriangular)
   {
      nNug = details::GeostrophicTools::cylTruncNug(nL, this->mIsTriangular);
   }
   else
   {
      nNug = details::GeostrophicTools::cylTruncNugC(nN, nL);
   }
   const auto nli = details::GeostrophicTools::nlist(nNug - 1, nL);
   const auto& nCpu = this->mNcpu;
   const auto& alpha = this->mUgAlpha;
   const auto& beta = this->mUgBeta;

   Internal::Array igridz;
   Internal::Array iweightz;
   details::GeostrophicTools::computeQuadratureZ(igridz, iweightz, nZ);

   if (alpha != 0.5_mp || beta != 1.0_mp)
   {
      throw std::logic_error("Generic (alpha,beta) pair is not implemented!");
   }

   // compute Gauss-Jacobi quadrature in x
   Internal::Array igridx, ilambda;
   Polynomial::Quadrature::JacobiRule jRule(alpha, beta - 1);
   jRule.computeQuadrature(igridx, ilambda, nS);

   // grid in s
   Internal::Array igrids =
      ((igridx.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

   // Compute integrator for geostrophic basis (weighted projector)
   Internal::Matrix iweights_proj;
   iweights_proj.resize(nS, nNug);
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
   wnl.compute<Internal::MHDFloat>(iweights_proj, nNug, 1, igrids,
      ilambda * c, Polynomial::Worland::Evaluator::Set());

   // allocate integration matrices
   mat.resize(nNug, nN * nL);
   mat.setConstant(0.0);

   // compute integration matrices
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
         details::GeostrophicTools::integrateZ(l, maxN_, iintgz, igrids, igridz,
            iweightz);
         Internal::Matrix tmp = iweights_proj.transpose() * iintgz;
         mat.block(0, l * nN, nNug, maxN_ + 1) = Internal::cast(tmp);
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
