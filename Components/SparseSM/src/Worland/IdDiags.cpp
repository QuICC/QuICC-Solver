/** 
 * @file IdDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/IdDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   IdDiags::IdDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l), mQ(q)
   {
   }

   IdDiags::ACoeff_t IdDiags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      if(this->mQ > 0)
      {
         val.topRows(this->mQ).setZero();
      }
      else if(this->mQ < 0)
      {
         val.bottomRows(-this->mQ).setZero();
      }

      return val;
   }

}
}
}
