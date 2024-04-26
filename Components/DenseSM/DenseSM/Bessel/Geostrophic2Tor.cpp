/**
 * @file Geostrophic2Tor.cpp
 * @brief Source of the implementation of geostrophic to toroidal basis conversion
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "Geostrophic2Tor.hpp"
#include "DenseSM/Bessel/details/GeostrophicTools.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Bessel/Generic.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   Geostrophic2Tor::Geostrophic2Tor(const int nN, const int nL, const int nCpu, const Scalar_t sDNu, const Scalar_t torDNu)
      : IMatrixSMOperator(nL*nN, nN), mNn(nN), mNl(nL), mNcpu(nCpu), mSDNu(sDNu), mTorDNu(torDNu)
   {
   }

   void Geostrophic2Tor::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      const auto& nN = this->mNn;
      const auto& nL = this->mNl;
      const auto& nR = DenseSM::Bessel::details::GeostrophicTools::cylTruncNr(nN, nL);
      int nNug = details::GeostrophicTools::cylTruncNug(nN, nL);

      // Create index list
      std::vector<int> nIdx;
      for (int n = 0; n < nNug; n++)
      {
         if (QuICCEnv().id() == n % this->mNcpu)
         {
            nIdx.push_back(n);
         }
      }

      // Create radial grid
      Internal::Array igrid, iweights;
      Polynomial::Quadrature::WorlandSphEnergyRule wquad;
      wquad.computeQuadrature(igrid, iweights, nR);

      //
      // Build operator
      mat = Internal::Matrix::Zero(nL*nN, nNug);

      std::map<int,Internal::Matrix> matP;
      for(int l = 0; l < nL; l++)
      {
         if(l % 2 == 1)
         {
            int nN_ = nN-1;
            Internal::Matrix ipoly(igrid.size(), nN_ + 1);
            Polynomial::Bessel::Generic<Polynomial::Bessel::SphJnl> jnl(this->mTorDNu);
            jnl.compute<Internal::MHDFloat>(ipoly, nN_ + 1, l, igrid, iweights);
            matP.try_emplace(l,ipoly);
         }
      }

      // Create Legendre grid and weights
      int nAlPoly = nL;
      int nTh = 3*(nAlPoly+1)/2;
      Internal::Array ialgrid, ialweights;
      Polynomial::Quadrature::LegendreRule lquad;
      lquad.computeQuadrature(ialgrid, ialweights, nTh);

      // Compute Legendre operator
      Polynomial::ALegendre::dPlm dplm;
      namespace evAL = Polynomial::ALegendre::Evaluator;
      Internal::Matrix alOp(ialgrid.size(), nAlPoly);
      dplm.compute<Internal::MHDFloat>(alOp, nAlPoly, 0, ialgrid, ialweights, evAL::Set());
      Internal::Array invLaplh = Internal::Array::LinSpaced(nAlPoly, 0, nAlPoly-1);
      invLaplh = (invLaplh.array()*(invLaplh.array() + 1.0)).pow(-1);
      invLaplh(0) = 0.0;

      using namespace Internal::Literals;
      // Create basis for geostrophic flow (used with l = 1)
      const auto& sDNu = this->mSDNu;
      Polynomial::Bessel::Generic<Polynomial::Bessel::SphJnl> ugJnl(sDNu);

      Internal::Array iugGrid(ialgrid.size());
      for(int n_: nIdx)
      {
         // Normalization (includes 2\pi from Fourier)
         Internal::MHDFloat scale = 2_mp*Internal::Math::PI;

         // Convert geostrophic flow into 2D spherical flow (r, l)
         // (spectral theta, ignore phi direction)
         Internal::Matrix tPoly = Internal::Matrix::Zero(igrid.size(), nAlPoly);
         for(int tk = 0; tk < igrid.size(); tk++)
         {
            // Convert cylidrical s to (r,theta)
            iugGrid = igrid(tk)*ialgrid.array().acos().sin();

            // Evaluate geostrophic flow in (r,theta)
            Internal::Matrix ipoly(iugGrid.size(), n_+1);
            ugJnl.compute<Internal::MHDFloat>(ipoly, n_+1, 1, iugGrid, Internal::Array());

            tPoly.row(tk) = -(invLaplh.asDiagonal()*(alOp.transpose()*(scale*ipoly.rightCols(1)))).transpose();
         }
         for(int l = 0; l < std::min(2*n_+2, nL); l++)
         {
            if(l%2 == 1)
            {
               // Compute spectral expansion
               Internal::Matrix tmp = matP.at(l).transpose()*(tPoly.col(l));

               mat.block(l*nN, n_, tmp.rows(), 1) = tmp;
            }
         }
      }

#if defined QUICC_MPI
      MPI_Allreduce(MPI_IN_PLACE, mat.data(), mat.size(),
         MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
   }

} // Bessel
} // DenseSM
} // QuICC
