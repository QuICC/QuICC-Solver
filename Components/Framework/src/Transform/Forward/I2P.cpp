/**
 * @file I2P.cpp
 * @brief Source of forward projection operator I2P
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2P.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2P::sTag()
   {
      return "Fwd::I2P";
   }

   const std::size_t& I2P::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2P>(I2P::sTag());
      return *i;
   }

   I2P::I2P()
      : IOperator(I2P::sTag())
   {
   }

   I2P::~I2P()
   {
   }

}
}
}
