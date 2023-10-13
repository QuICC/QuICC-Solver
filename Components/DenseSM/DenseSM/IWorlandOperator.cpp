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
#include "DenseSM/IWorlandOperator.hpp"
#include "DenseSM/Worland/Tools.hpp"

namespace QuICC {

namespace DenseSM {

   IWorlandOperator::IWorlandOperator(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta)
      : IDenseSMOperator(rows, cols)
   {
      this->mType = Worland::Tools::identifyBasis(alpha, dBeta);
   }

   Worland::WorlandKind IWorlandOperator::type() const
   {
      return this->mType;
   }

}
}
