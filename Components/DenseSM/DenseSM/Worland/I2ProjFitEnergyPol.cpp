/**
 * @file I2ProjFitEnergyPol.cpp
 * @brief Source of the implementation of the full sphere Worland projection operator onto best energy fit for poloidal scalar with I2 quasi-inverse
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "I2ProjFitEnergyPol.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   I2ProjFitEnergyPol::I2ProjFitEnergyPol(const int outRows, const std::size_t bcId, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta), mOutRows(outRows), mL(l), mProj(outRows, bcId, rows, cols, alpha, dBeta, l, q)
   {
   }

   void I2ProjFitEnergyPol::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
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

   void I2ProjFitEnergyPol::buildChebyshevOp(Internal::Matrix& mat, const int rows, const int cols) const
   {
      const auto& nbar = this->mOutRows;

      if(nbar > 1)
      {
         const auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
         const auto db = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
         const auto& l = this->mL;

         Polynomial::Quadrature::WorlandLegendreRule wquad;

         mat = this->mProj.mat();

         SparseSM::Worland::I2 i2(mat.rows(), mat.rows(), a, db, l);
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
