/**
 * @file I2ProjFitEnergyPol.cpp
 * @brief Source of the implementation of the full sphere Worland projection operator onto best energy fit for poloidal scalar with I2 quasi-inverse
 */

// System includes
//

// Project includes
//
#include "I2ProjFitEnergyPol.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   I2ProjFitEnergyPol::I2ProjFitEnergyPol(const int outRows, const std::size_t bcId, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IEmbeddedOperator(rows, cols, alpha, dBeta), mOutRows(outRows), mL(l), mProj(outRows, bcId, rows, cols, alpha, dBeta, l, q)
   {
   }

   void I2ProjFitEnergyPol::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
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

   void I2ProjFitEnergyPol::buildGenericOp(Internal::Matrix& mat, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta) const
   {
      const auto& nbar = this->mOutRows;

      if(nbar > 1)
      {
         const auto& l = this->mL;

         mat = this->mProj.mat();

         SparseSM::Worland::I2 i2(mat.rows(), mat.rows(), alpha, dBeta, l);
         mat = i2.mat()*mat;
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
