/**
 * @file ProjFitEnergyPol.cpp
 * @brief Source of the implementation of the full sphere Worland projection operator onto best energy fit for poloidal scalar
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "ProjFitEnergyPol.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphere.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   ProjFitEnergyPol::ProjFitEnergyPol(const int outRows, const std::size_t bcId, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta), mOutRows(outRows), mL(l), mBcId(bcId)
   {
   }

   void ProjFitEnergyPol::buildOpImpl(internal::Matrix& mat, const int rows, const int cols) const
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

   void ProjFitEnergyPol::buildChebyshevOp(internal::Matrix& mat, const int rows, const int cols) const
   {
      const auto& nbar = this->mOutRows;

      if(nbar > 1)
      {
         const auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
         const auto db = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
         const auto& l = this->mL;
         const auto& dl = static_cast<internal::MHDFloat>(this->mL);
         const auto ll1 = dl*(dl+MHD_MP(1.0));
         Polynomial::Quadrature::WorlandLegendreRule wquad;

         int rp = 2*rows + l;
         int pts = rp + 2 + (rp + 2)%2;
         internal::Array igrid;
         internal::Array iweights;
         wquad.computeQuadrature(igrid, iweights, pts);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         Polynomial::Worland::drWnl drwnl;

         internal::Matrix tmpBwd(igrid.size(), rows);
         wnl.compute<internal::MHDFloat>(tmpBwd, rows, l, igrid, internal::Array(), ev::Set());
         internal::Matrix tmpFwd(igrid.size(), rows);
         wnl.compute<internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array(), ev::Set());
         internal::Matrix matW = ll1*ll1*(tmpFwd.transpose()*tmpBwd);

         drwnl.compute<internal::MHDFloat>(tmpBwd, rows, l, igrid, internal::Array(), ev::Set());
         drwnl.compute<internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array(), ev::Set());
         matW += ll1*(tmpFwd.transpose()*tmpBwd);

         SparseSM::Worland::Boundary::Value bc(a, db, l);
         internal::Matrix bcVal = bc.compute(rows-1).matrix();
         internal::Matrix bcMat = (dl*ll1*bcVal)*bcVal.transpose();

         matW += bcMat;

         SparseMatrix matS;
         if(this->mBcId == Bc::Name::Insulating::id())
         {
            SparseSM::Worland::Stencil::InsulatingSphere S(nbar, nbar-1, a, db, l);
            matS = S.mat();
         }
         else
         {
            throw std::logic_error("Unknown boundary condition");
         }

         internal::Matrix matWbar = (matS.transpose()*matW.block(0,0,nbar,nbar)*matS);
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
