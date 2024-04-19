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
#include "DenseSM/Bessel/details/GeostrophicTools.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Tor2GridS.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

Tor2GridS::Tor2GridS(const int nN, const int nL, const int nCpu, const Scalar_t sDNu,
   const Scalar_t torDNu) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mNcpu(nCpu),
    mSDNu(sDNu),
    mTorDNu(torDNu)
{}

void Tor2GridS::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
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
   iweightS.resize(0);

   // allocate integration matrices
   mat.resize(nS, nN * nL);
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
         details::GeostrophicTools::integrateZ(l, maxN_, iintgz, igridS, igridZ, iweightZ, this->mTorDNu);
         mat.block(0, l * nN, nS, maxN_ + 1) = iintgz;
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
