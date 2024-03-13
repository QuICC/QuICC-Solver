/**
 * @file ProjFitEnergyR2.cpp
 * @brief Source of the implementation of the full sphere Worland projection operator onto best energy fit
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "ProjFitEnergyR2.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   ProjFitEnergyR2::ProjFitEnergyR2(const int outRows, const std::size_t bcId, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IEmbeddedOperator(rows, cols, alpha, dBeta), mOutRows(outRows), mL(l), mBcId(bcId)
   {
   }

   void ProjFitEnergyR2::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
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

   void ProjFitEnergyR2::buildChebyshevOp(Internal::Matrix& mat, const int rows, const int cols) const
   {
      const auto& nbar = this->mOutRows;

      if(nbar > 1)
      {
         const Polynomial::Worland::worland_chebyshev_t wt;
         const auto a = wt.ALPHA;
         const auto db =wt.DBETA;
         const auto& l = this->mL;
         Polynomial::Quadrature::WorlandLegendreRule wquad;

         int rp = 2*rows + l;
         int pts = rp + 2 + (rp + 2)%2;
         Internal::Array igrid;
         Internal::Array iweights;
         wquad.computeQuadrature(igrid, iweights, pts);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         Internal::Matrix tmpBwd(igrid.size(), rows);
         wnl.compute<Internal::MHDFloat>(tmpBwd, rows, l, igrid, Internal::Array(), ev::Set());
         Internal::Matrix tmpFwd(igrid.size(), rows);
         wnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array()*igrid.array().abs2(), ev::Set());
         Internal::Matrix matW = tmpFwd.transpose()*tmpBwd;

         SparseMatrix matS;
         if(this->mBcId == Bc::Name::FixedTemperature::id() || this->mBcId == Bc::Name::Insulating::id())
         {
            SparseSM::Worland::Stencil::Value S(nbar, nbar-1, a, db, l);
            matS = S.mat();
         }
         else if(this->mBcId == Bc::Name::FixedFlux::id())
         {
            SparseSM::Worland::Stencil::D1 S(nbar, nbar-1, a, db, l);
            matS = S.mat();
         }
         else
         {
            throw std::logic_error("Unknown boundary condition");
         }


         Internal::Matrix matWbar = (matS.transpose()*matW.block(0,0,nbar,nbar)*matS);
         mat.resize(rows, cols);
         mat.topRows(nbar) = matS*matWbar.inverse()*matS.transpose()*matW.topRows(nbar);
         mat.bottomRows(rows-nbar).setZero();
      }
      else
      {
         mat.setZero(rows, rows);
      }
   }

} // Worland
} // DenseSM
} // QuICC
