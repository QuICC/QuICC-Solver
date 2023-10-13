/**
 * @file IWorlandOperator.cpp
 * @brief Source of the implementation of generic interface to a full sphere Worland sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <unsupported/Eigen/SpecialFunctions>

// Project includes
//
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/Tools.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace SparseSM {

   IWorlandOperator::IWorlandOperator(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta)
      : ISparseSMOperator(rows, cols)
   {
      this->mType = Worland::Tools::identifyBasis(alpha, dBeta);
   }

   Worland::WorlandKind IWorlandOperator::type() const
   {
      return this->mType;
   }

}
}
