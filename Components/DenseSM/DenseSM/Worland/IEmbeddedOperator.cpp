/** 
 * @file IEmbeddedOperator.cpp
 * @brief Source of the implementation of generic interface to a full sphere Worland dense operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/IEmbeddedOperator.hpp"
#include "DenseSM/Worland/Tools.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   IEmbeddedOperator::IEmbeddedOperator(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta)
      : DenseSM::IEmbeddedSMOperator(rows, cols)
   {
      this->mType = Worland::Tools::identifyBasis(alpha, dBeta);
   }

   Worland::WorlandKind IEmbeddedOperator::type() const
   {
      return this->mType;
   }

}
}
}
