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
#include "Geostrophic2Tor.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandChebyshevRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   Geostrophic2Tor::Geostrophic2Tor(const int nN, const int maxnl, const int nR, const int maxNug, const ArrayI& nli, const std::vector<int>& nIdx, const Scalar_t ugAlpha, const Scalar_t ugDBeta, const Scalar_t alpha, const Scalar_t dBeta, const int q)
      : IGeostrophicOperator(ugAlpha, ugDBeta, maxnl*nN, maxNug+1, alpha, dBeta, q), mNn(nN), mMaxnl(maxnl), mNr(nR), mMaxNug(maxNug), mNlist(nli), mNidx(nIdx)
   {
   }

   void Geostrophic2Tor::buildOpImpl(internal::Matrix& mat, const int rows, const int cols) const
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->buildChebyshevOp(mat, rows, cols);
            break;
         case WorlandKind::LEGENDRE:
            throw std::logic_error("Legendre basis operator not implemented");
            break;
         case WorlandKind::CYLENERGY:
            throw std::logic_error("Cylindrical energy basis operator not implemented");
            break;
         case WorlandKind::SPHENERGY:
            throw std::logic_error("Spherical energy basis operator not implemented");
            break;
      }
   }

   void Geostrophic2Tor::buildChebyshevOp(internal::Matrix& mat, const int rows, const int cols) const
   {
      const auto& nN = this->mNn;
      const auto& maxnl = this->mMaxnl;
      const auto& nR = this->mNr;
      const auto& maxNug = this->mMaxNug;
      const auto& nli = this->mNlist;
      const auto& nIdx = this->mNidx;
      const auto& a = this->mcUgAlpha;
      const auto& b = this->mcUgDBeta;

      internal::Array igrid, iweights;
      Polynomial::Quadrature::WorlandChebyshevRule wquad;
      wquad.computeQuadrature(igrid, iweights, nR);

      mat = internal::Matrix::Zero(maxnl*nN, maxNug+1);

      std::map<int,internal::Matrix> matP;
      for(int l = 0; l < maxnl; l++)
      {
         if(l % 2 == 1)
         {
            int nN_ = std::min(nli(l), nN-1);
            internal::Matrix ipoly(igrid.size(), nN_ + 1);
            Polynomial::Worland::Wnl wnl;
            wnl.compute<internal::MHDFloat>(ipoly, nN_ + 1, l, igrid, iweights, Polynomial::Worland::Evaluator::Set());
            matP.try_emplace(l,ipoly);
         }
      }

      // Create Legendre grid and weights
      int nAlPoly = maxnl;
      int nTh = 3*(nAlPoly+1)/2;
      internal::Array ialgrid, ialweights;
      Polynomial::Quadrature::LegendreRule lquad;
      lquad.computeQuadrature(ialgrid, ialweights, nTh);

      // Compute Legendre operator
      Polynomial::ALegendre::dPlm dplm;
      namespace evAL = Polynomial::ALegendre::Evaluator;
      internal::Matrix alOp(ialgrid.size(), nAlPoly);;
      dplm.compute<internal::MHDFloat>(alOp, nAlPoly, 0, ialgrid, ialweights, evAL::Set());
      internal::Array invLaplh = internal::Array::LinSpaced(nAlPoly, 0, nAlPoly-1);
      invLaplh = (invLaplh.array()*(invLaplh.array() + 1.0)).pow(-1);
      invLaplh(0) = 0.0;

      // Create basis for geostrophic flow (used with l = 1)
      internal::MHDFloat ugA;
      internal::MHDFloat ugDB;
      if(this->isUgBasis(this->mcUgAlpha, this->mcUgDBeta))
      {
         // Using \tilde{\Lambda}(s) basis
         ugA = this->mcUgAlpha;
         ugDB = this->mcUgDBeta;
         throw std::logic_error("General geostrophic basis is not implemented");
      }
      else
      {
         // Using \Lambda(s) basis
         ugA = MHD_MP(0.5);
         ugDB = MHD_MP(0);
      }
      internal::MHDFloat ugB = ugDB + MHD_MP(1);
      Polynomial::Worland::Wnl ugWnl(ugA, ugDB);

      internal::Array iugGrid(ialgrid.size());
      for(int n_: nIdx)
      {
         // Normalization (includes 2\pi from Fourier)
         internal::MHDFloat scale = this->Cn(n_)*(MHD_MP(2)*Precision::PI)/this->Cnab(n_, ugA, ugB);

         // Convert geostrophic flow into 2D spherical flow (r, l)
         // (spectral theta, ignore phi direction)
         internal::Matrix tPoly = internal::Matrix::Zero(igrid.size(), nAlPoly);
         for(int tk = 0; tk < igrid.size(); tk++)
         {
            // Convert cylidrical s to (r,theta)
            iugGrid = igrid(tk)*ialgrid.array().acos().sin();

            // Evaluate geostrophic flow in (r,theta)
            internal::Matrix ipoly(iugGrid.size(), n_+1);
            ugWnl.compute<internal::MHDFloat>(ipoly, n_+1, 1, iugGrid, internal::Array(), Polynomial::Worland::Evaluator::Set());

            tPoly.row(tk) = -(invLaplh.asDiagonal()*(alOp.transpose()*(scale*ipoly.rightCols(1)))).transpose();
         }
         for(int l = 0; l < std::min(2*n_+2, maxnl); l++)
         {
            if(l%2 == 1)
            {
               // Compute Worland expansion
               internal::Matrix tmp = matP.at(l).transpose()*(tPoly.col(l));

               mat.block(l*nN, n_, tmp.rows(), 1) = tmp;
            }
         }
      }
   }

} // Worland
} // DenseSM
} // QuICC
