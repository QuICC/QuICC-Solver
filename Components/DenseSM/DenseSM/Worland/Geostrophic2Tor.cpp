/**
 * @file Geostrophic2Tor.cpp
 * @brief Source of the implementation of the full sphere Worland projection operator onto best energy fit
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
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   Geostrophic2Tor::Geostrophic2Tor(const int nN, const int nL, const int nCpu, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const bool isTriangular)
      : IGeostrophicOperator(ugAlpha, ugBeta, isGenericBasis, nL*nN, nN), mNn(nN), mNl(nL), mNcpu(nCpu), mIsTriangular(isTriangular)
   {
   }

   Geostrophic2Tor::Geostrophic2Tor(const int nN, const int nL, const int nCpu, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const bool isTriangular, const Scalar_t alpha, const Scalar_t dBeta)
      : IGeostrophicOperator(ugAlpha, ugBeta, isGenericBasis, nL*nN, nN, alpha, dBeta), mNn(nN), mNl(nL), mNcpu(nCpu), mIsTriangular(isTriangular)
   {
   }

   void Geostrophic2Tor::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      const auto& nN = this->mNn;
      const auto& nL = this->mNl;
      const auto& nR = DenseSM::Worland::details::GeostrophicTools::cylTruncNr(nL, this->mIsTriangular);
      int nNug;
      if(this->mIsTriangular)
      {
         nNug = details::GeostrophicTools::cylTruncNug(nL, this->mIsTriangular);
      }
      else
      {
         nNug = details::GeostrophicTools::cylTruncNugC(nN, nL);
      }

      // Create truncation list
      const auto nli = details::GeostrophicTools::nlist(nNug - 1, nL);

      // Create index list
      std::vector<int> nIdx;
      for (int n = 0; n < nNug; n++)
      {
         if (QuICCEnv().id() == n % this->mNcpu)
         {
            nIdx.push_back(n);
         }
      }

      // Select Worland type and create quadrature grid and weights
      Scalar_t alpha, dBeta;
      Internal::Array igrid, iweights;

      auto define_worland = [&](auto& w)
      {
         alpha = w.ALPHA;
         dBeta = w.DBETA;
         typename std::remove_reference<decltype(w)>::type::Rule quad;
         quad.computeQuadrature(igrid, iweights, nR);
      };

      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            {
               ::QuICC::Polynomial::Worland::worland_chebyshev_t wt;
               define_worland(wt);
               break;
            }
         case WorlandKind::LEGENDRE:
            {
               ::QuICC::Polynomial::Worland::worland_legendre_t wt;
               define_worland(wt);
               break;
            }
         case WorlandKind::CYLENERGY:
            {
               ::QuICC::Polynomial::Worland::worland_cylenergy_t wt;
               define_worland(wt);
               break;
            }
         case WorlandKind::SPHENERGY:
            {
               ::QuICC::Polynomial::Worland::worland_sphenergy_t wt;
               define_worland(wt);
               break;
            }
      }

      //
      // Build operator
      mat = Internal::Matrix::Zero(nL*nN, nNug);

      std::map<int,Internal::Matrix> matP;
      for(int l = 0; l < nL; l++)
      {
         if(l % 2 == 1)
         {
            int nN_ = std::min(nli(l), nN-1);
            Internal::Matrix ipoly(igrid.size(), nN_ + 1);
            Polynomial::Worland::Wnl wnl(alpha, dBeta);
            wnl.compute<Internal::MHDFloat>(ipoly, nN_ + 1, l, igrid, iweights, Polynomial::Worland::Evaluator::Set());
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
      const auto& ugA = this->mcUgAlpha;
      const auto& ugB = this->mcUgBeta;
      const auto ugDB = ugB - 1_mp;
      if(this->mcIsGenericBasis)
      {
         throw std::logic_error("General geostrophic basis is not implemented");
      }
      Polynomial::Worland::Wnl ugWnl(ugA, ugDB);

      Internal::Array iugGrid(ialgrid.size());
      for(int n_: nIdx)
      {
         // Normalization (includes 2\pi from Fourier)
         Internal::MHDFloat scale = details::GeostrophicTools::Cn(n_)*(2_mp*Internal::Math::PI)/details::GeostrophicTools::Cnab(n_, ugA, ugB);

         // Convert geostrophic flow into 2D spherical flow (r, l)
         // (spectral theta, ignore phi direction)
         Internal::Matrix tPoly = Internal::Matrix::Zero(igrid.size(), nAlPoly);
         for(int tk = 0; tk < igrid.size(); tk++)
         {
            // Convert cylidrical s to (r,theta)
            iugGrid = igrid(tk)*ialgrid.array().acos().sin();

            // Evaluate geostrophic flow in (r,theta)
            Internal::Matrix ipoly(iugGrid.size(), n_+1);
            ugWnl.compute<Internal::MHDFloat>(ipoly, n_+1, 1, iugGrid, Internal::Array(), Polynomial::Worland::Evaluator::Set());

            tPoly.row(tk) = -(invLaplh.asDiagonal()*(alOp.transpose()*(scale*ipoly.rightCols(1)))).transpose();
         }
         for(int l = 0; l < std::min(2*n_+2, nL); l++)
         {
            if(l%2 == 1)
            {
               // Compute Worland expansion
               Internal::Matrix tmp = matP.at(l).transpose()*(tPoly.col(l));

               mat.block(l*nN, n_, tmp.rows(), 1) = tmp;
            }
         }
      }

#if defined QUICC_MPI
      MPI_Allreduce(MPI_IN_PLACE, this->mGeo2Tor.data(), this->mGeo2Tor.size(),
         MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
   }

} // Worland
} // DenseSM
} // QuICC
