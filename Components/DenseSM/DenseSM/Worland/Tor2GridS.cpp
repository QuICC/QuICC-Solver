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
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Tor2GridS.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Casts.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

Tor2GridS::Tor2GridS(const int nN, const int nL, const int nCpu, const Scalar_t alpha,
   const Scalar_t beta, const bool isGenericBasis, const Scalar_t alphaB,
   const Scalar_t betaB, const bool isTriangular) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mNcpu(nCpu),
    mUgAlpha(alpha),
    mUgBeta(beta),
    mIsGenericBasis(isGenericBasis),
    mAlphaB(alphaB),
    mBetaB(betaB),
    mIsTriangular(isTriangular)
{}

void Tor2GridS::buildOpImpl(Internal::Matrix& mat, const int rows,
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
   const auto nli = details::GeostrophicTools::nlist(nNug-1, nL);
   const auto& nCpu = this->mNcpu;
   const auto& alpha = this->mUgAlpha;
   const auto& beta = this->mUgBeta;
   const auto& alphaB = this->mAlphaB;
   const auto& betaB = this->mBetaB;

   Internal::Array igridz;
   Internal::Array iweightz;
   details::GeostrophicTools::computeQuadratureZ(igridz, iweightz, nZ);

   if (alpha != 0.5_mp || beta != 1.0_mp)
   {
      throw std::logic_error("Generic (alpha,beta) pair is not implemented!");
   }

   // compute Gauss-Jacobi quadrature in x
   Internal::Array igridx, ilambda;
   Polynomial::Quadrature::JacobiRule jRuleA(alpha, beta - 1);
   jRuleA.computeQuadrature(igridx, ilambda, nS);

   // grid in s
   Internal::Array igrids =
      ((igridx.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

   // compute Gauss-Jacobi quadrature in x for second alpha,beta pair
   Polynomial::Quadrature::JacobiRule jRuleB(alphaB, betaB);
   jRuleB.computeQuadrature(igridx, ilambda, nS);

   // grid in s
   Internal::Array igridsB =
      ((igridx.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

   // allocate integration matrices
   mat.resize(nS, nN * nL);
   mat.setConstant(0.0);

   // compute integration matrices
   int pid = 0;
   for (int l = 1; l < nL; l += 2)
   {
      if (QuICCEnv().id() == (pid % nCpu))
      {
         int maxN_;
         Internal::Array* pGrids;
         if (this->mIsTriangular)
         {
            maxN_ = nli(l);
            if (this->mIsGenericBasis)
            {
               pGrids = &igridsB;
            }
            else
            {
               pGrids = &igrids;
            }
         }
         else
         {
            maxN_ = nN - 1;
            pGrids = &igrids;
         }

         Internal::Matrix iintgz;
         details::GeostrophicTools::integrateZ(l, maxN_, iintgz,                  *pGrids, igridz, iweightz);
         mat.block(0, l * nN, nS, maxN_ + 1) = Internal::cast(iintgz);
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
