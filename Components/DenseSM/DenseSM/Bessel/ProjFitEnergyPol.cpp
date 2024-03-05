/**
 * @file ProjFitEnergyPol.cpp
 * @brief Source of the implementation of the full sphere Bessel projection operator onto best energy fit for poloidal scalar
 */

// System includes
//
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "ProjFitEnergyPol.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"
#include "QuICC/Polynomial/Bessel/Insulating.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Polynomial/Bessel/drSphJnl.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   ProjFitEnergyPol::ProjFitEnergyPol(const int outRows, const std::size_t bcId, const int rows, const int cols, const int l)
      : IEmbeddedSMOperator(rows, cols), mOutRows(outRows), mL(l), mBcId(bcId)
   {
   }

   void ProjFitEnergyPol::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      this->buildGenericOp(mat, rows, cols);
   }

   void ProjFitEnergyPol::buildGenericOp(Internal::Matrix& mat, const int rows, const int cols) const
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

         Polynomial::Bessel::Insulating<Polynomial::Bessel::SphJnl> jnl;
         Polynomial::Bessel::Insulating<Polynomial::Bessel::drSphJnl> drjnl;

         Internal::Matrix tmpBwd(igrid.size(), rows);
         jnl.compute<Internal::MHDFloat>(tmpBwd, rows, l, igrid, Internal::Array());
         Internal::Matrix tmpFwd(igrid.size(), rows);
         jnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array());
         Internal::Matrix matW = ll1*ll1*(tmpFwd.transpose()*tmpBwd);

         drjnl.compute<Internal::MHDFloat>(tmpBwd, rows, l, igrid, Internal::Array());
         drjnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array());
         matW += ll1*(tmpFwd.transpose()*tmpBwd);

         Internal::Array iend = Internal::Array::Ones(1);
         Internal::Matrix bcVal(1, rows);
         jnl.compute<Internal::MHDFloat>(bcVal, rows, l, iend, Internal::Array());
         Internal::Matrix bcMat = (dl*ll1*bcVal.transpose())*bcVal;

         matW += bcMat;

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
