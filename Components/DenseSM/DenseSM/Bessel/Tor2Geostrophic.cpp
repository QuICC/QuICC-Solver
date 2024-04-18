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
#include "DenseSM/Bessel/details/GeostrophicTools.hpp"
#include "QuICC/Polynomial/Bessel/Generic.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

Tor2Geostrophic::Tor2Geostrophic(const int nN, const int nL, const int nCpu,
   const Scalar_t sDNu, const Scalar_t torDNu) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mNcpu(nCpu),
    mSDNu(sDNu),
    mTorDNu(torDNu)
{}

void Tor2Geostrophic::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   using namespace Internal::Literals;
   const auto& nN = this->mNn;
   const auto& nL = this->mNl;

   // Compute Z integral quadrature
   const auto nZ = details::GeostrophicTools::cylTruncNz(nN, nL);
   Internal::Array igridZ;
   Internal::Array iweightZ;
   details::GeostrophicTools::computeQuadratureZ(igridZ, iweightZ, nZ);

   // Compute S grid
   const auto nS = details::GeostrophicTools::cylTruncNs(nN, nL);
   const auto& sDNu = this->mSDNu;
   Internal::Array igridS;
   Internal::Array iweightS;
   details::GeostrophicTools::computeGridS(igridS, iweightS, nS, sDNu);

   // Compute integrator for geostrophic basis (weighted projector)
   int nNug = details::GeostrophicTools::cylTruncNug(nN, nL);
   Internal::Matrix iweightS_proj;
   iweightS_proj.resize(nS, nNug);
   Polynomial::Bessel::Generic<Polynomial::Bessel::SphJnl> jnl(sDNu);
   jnl.compute<Internal::MHDFloat>(iweightS_proj, nNug, 1, igridS,
      iweightS);

   // allocate integration matrices
   mat.resize(nNug, nN * nL);
   mat.setConstant(0.0);

   // compute integration matrices
   const auto& nCpu = this->mNcpu;
   int pid = 0;
   for (int l = 1; l < nL; l += 2)
   {
      if (QuICCEnv().id() == (pid % nCpu))
      {
         int maxN_ = nN - 1;

         Internal::Matrix iintgz;
         details::GeostrophicTools::integrateZ(l, maxN_, iintgz, igridS, igridZ,
            iweightZ, this->mTorDNu);
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

} // namespace Bessel
} // namespace DenseSM
} // namespace QuICC
