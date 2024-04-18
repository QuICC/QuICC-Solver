/**
 * @file Tor2GridS.cpp
 * @brief Source of the implementation of the projection operator from the
 * toroidal scalar to the geostrophic basis
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Tor2GridS.hpp"
#include "Types/Internal/Casts.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

Tor2GridS::Tor2GridS(const int nN, const int nL, const int nCpu, const Scalar_t alpha,
   const Scalar_t beta, const bool isTriangular, const Scalar_t wAlpha, const Scalar_t wDBeta) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mNcpu(nCpu),
    mUgAlpha(alpha),
    mUgBeta(beta),
    mIsTriangular(isTriangular),
    mAlpha(wAlpha),
    mDBeta(wDBeta)
{}

Tor2GridS::Tor2GridS(const int nN, const int nL, const int nCpu, const Scalar_t alpha,
   const Scalar_t beta, const bool isTriangular) :
    Tor2GridS(nN, nL, nCpu, alpha, beta, isTriangular, Polynomial::Worland::worland_default_t::ALPHA, Polynomial::Worland::worland_default_t::DBETA)
{}

void Tor2GridS::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
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
   iweightS.resize(0);

   // allocate integration matrices
   mat.resize(nS, nN * nL);
   mat.setConstant(0.0);

   // compute integration matrices
   int nNug;
   if(this->mIsTriangular)
   {
      nNug = details::GeostrophicTools::cylTruncNug(nL, this->mIsTriangular);
   }
   else
   {
      nNug = details::GeostrophicTools::cylTruncNugC(nN, nL);
   }
   const auto nli = details::GeostrophicTools::nlist(nNug-1, nL);
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
         details::GeostrophicTools::integrateZ(l, maxN_, iintgz, igridS, igridZ, iweightZ, this->mAlpha, this->mDBeta);
         mat.block(0, l * nN, nS, maxN_ + 1) = iintgz;
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
