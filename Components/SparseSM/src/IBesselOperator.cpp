/**
 * @file IBesselOperator.cpp
 * @brief Source of the implementation of generic interface to a full sphere Bessel sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/IBesselOperator.hpp"

namespace QuICC {

namespace SparseSM {

   IBesselOperator::IBesselOperator(const int rows, const int cols, const Bessel::BesselKind type)
      : ISparseSMOperator(rows, cols), mType(type)
   {
   }

   Bessel::BesselKind IBesselOperator::type() const
   {
      return this->mType;
   }

}
}
