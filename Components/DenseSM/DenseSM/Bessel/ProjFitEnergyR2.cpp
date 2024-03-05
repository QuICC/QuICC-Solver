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
#include "QuICC/Polynomial/Bessel/Value.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   ProjFitEnergyR2::ProjFitEnergyR2(const int outRows, const std::size_t bcId, const int rows, const int cols, const int l)
      : IEmbeddedSMOperator(rows, cols), mOutRows(outRows), mL(l), mBcId(bcId)
   {
   }

   void ProjFitEnergyR2::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      this->buildGenericOp(mat, rows, cols);
   }

   void ProjFitEnergyR2::buildGenericOp(Internal::Matrix& mat, const int rows, const int cols) const
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

         Polynomial::Bessel::Value<Polynomial::Bessel::SphJnl> jnl;
         Internal::Matrix tmpBwd(igrid.size(), rows);
         jnl.compute<Internal::MHDFloat>(tmpBwd, rows, l, igrid, Internal::Array());
         Internal::Matrix tmpFwd(igrid.size(), rows);
         jnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array()*igrid.array().abs2());
         Internal::Matrix matW = tmpFwd.transpose()*tmpBwd;

         Internal::Matrix matWbar = (matW.block(0,0,nbar,nbar));
         mat.resize(rows, cols);
         mat.topRows(nbar) = matWbar.inverse()*matW.topRows(nbar);
         mat.bottomRows(rows-nbar).setZero();
      }
      else
      {
         mat.setZero(rows, rows);
      }
   }

} // namespace Bessel
} // namespace DenseSM
} // namespace QuICC
