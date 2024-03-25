/** 
 * @file IWorlandOperator.cpp
 * @brief Source of the implementation of generic interface to a full sphere Worland dense operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "DenseSM/Worland/Tools.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   IWorlandOperator::IWorlandOperator(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta)
      : DenseSM::IMatrixSMOperator(rows, cols), mcAlpha(alpha), mcDBeta(dBeta)
   {
      this->mType = Worland::Tools::identifyBasis(this->mcAlpha, this->mcDBeta);
   }

   IWorlandOperator::IWorlandOperator(const int rows, const int cols)
      : IWorlandOperator(rows, cols, Polynomial::Worland::worland_default_t::ALPHA, Polynomial::Worland::worland_default_t::DBETA)
   {}

   Worland::WorlandKind IWorlandOperator::type() const
   {
      return this->mType;
   }

}
}
}
