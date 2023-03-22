/** 
 * @file ICondition.cpp
 * @brief Source of the base implementation for a boundary condition
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <limits>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   ICondition::ICondition(const Scalar_t lower, const Scalar_t upper, const Position pos)
      : mPosition(pos), mA(std::numeric_limits<Scalar_t>::max()), mB(std::numeric_limits<Scalar_t>::max())
   {
      this->setBounds(lower,upper);
   }

   void ICondition::setBounds(const Scalar_t lower, const Scalar_t upper)
   {
      this->mA = (upper - lower)/MHD_MP(2.0);
      this->mB = (upper + lower)/MHD_MP(2.0);
   }

   ICondition::Position ICondition::position() const
   {
      return this->mPosition;
   }

   ICondition::Scalar_t ICondition::a() const
   {
      return this->mA;
   }

   ICondition::Scalar_t ICondition::b() const
   {
      return this->mB;
   }

   ICondition::Scalar_t ICondition::c() const
   {
      return 2.0;
   }

} // Boundary
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
