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
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
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
      Scalar_t alpha, dBeta;
      auto define_worland = [&](auto& w)
      {
         alpha = w.ALPHA;
         dBeta = w.DBETA;
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
      this->buildGenericOp(mat, rows, cols, alpha, dBeta);
   }

   void ProjFitEnergyR2::buildGenericOp(Internal::Matrix& mat, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta) const
   {
      const auto& nbar = this->mOutRows;

      if(nbar > 1)
      {
         const auto& l = this->mL;

         // Compute Legendre quadrature to integrate r polynomial
         int rp = 2*rows + l;
         int pts = rp + 2 + (rp + 2)%2;
         Internal::Array igrid;
         Internal::Array iweights;
         Polynomial::Quadrature::WorlandLegendreRule wquad;
         wquad.computeQuadrature(igrid, iweights, pts);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl(alpha, dBeta);
         Internal::Matrix tmpBwd(igrid.size(), rows);
         wnl.compute<Internal::MHDFloat>(tmpBwd, rows, l, igrid, Internal::Array(), ev::Set());
         Internal::Matrix tmpFwd(igrid.size(), rows);
         wnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array()*igrid.array().abs2(), ev::Set());
         Internal::Matrix matW = tmpFwd.transpose()*tmpBwd;

         SparseMatrix matS;
         if(this->mBcId == Bc::Name::FixedTemperature::id() || this->mBcId == Bc::Name::Insulating::id())
         {
            SparseSM::Worland::Stencil::Value S(nbar, nbar-1, alpha, dBeta, l);
            matS = S.mat();
         }
         else if(this->mBcId == Bc::Name::FixedFlux::id())
         {
            SparseSM::Worland::Stencil::D1 S(nbar, nbar-1, alpha, dBeta, l);
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
