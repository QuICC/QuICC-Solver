/**
 * @file ProjFitEnergyPol.cpp
 * @brief Source of the implementation of the full sphere Worland projection operator onto best energy fit for poloidal scalar
 */

// System includes
//
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "ProjFitEnergyPol.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphere.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   ProjFitEnergyPol::ProjFitEnergyPol(const int outRows, const std::size_t bcId, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta), mOutRows(outRows), mL(l), mBcId(bcId)
   {
   }

   void ProjFitEnergyPol::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
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

   void ProjFitEnergyPol::buildGenericOp(Internal::Matrix& mat, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta) const
   {
      const auto& nbar = this->mOutRows;

      if(nbar > 1)
      {
         const auto& l = this->mL;
         const auto dl = static_cast<Internal::MHDFloat>(this->mL);
         const auto ll1 = dl*(dl+MHD_MP(1.0));

         // Compute Legendre quadrature to integrate r polynomial
         int rp = 2*rows + l;
         int pts = rp + 2 + (rp + 2)%2;
         Internal::Array igrid;
         Internal::Array iweights;
         Polynomial::Quadrature::WorlandLegendreRule wquad;
         wquad.computeQuadrature(igrid, iweights, pts);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl(alpha, dBeta);
         Polynomial::Worland::drWnl drwnl(alpha, dBeta);

         Internal::Matrix tmpBwd(igrid.size(), rows);
         wnl.compute<Internal::MHDFloat>(tmpBwd, rows, l, igrid, Internal::Array(), ev::Set());
         Internal::Matrix tmpFwd(igrid.size(), rows);
         wnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array(), ev::Set());
         Internal::Matrix matW = ll1*ll1*(tmpFwd.transpose()*tmpBwd);

         drwnl.compute<Internal::MHDFloat>(tmpBwd, rows, l, igrid, Internal::Array(), ev::Set());
         drwnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array(), ev::Set());
         matW += ll1*(tmpFwd.transpose()*tmpBwd);

         SparseSM::Worland::Boundary::Value bc(alpha, dBeta, l);
         Internal::Matrix bcVal = bc.compute(rows-1).matrix();
         Internal::Matrix bcMat = (dl*ll1*bcVal)*bcVal.transpose();

         matW += bcMat;

         SparseMatrix matS;
         if(this->mBcId == Bc::Name::Insulating::id())
         {
            SparseSM::Worland::Stencil::InsulatingSphere S(nbar, nbar-1, alpha, dBeta, l);
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
