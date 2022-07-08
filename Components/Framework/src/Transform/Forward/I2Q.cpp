/**
 * @file I2Q.cpp
 * @brief Source of forward projection operator I2Q
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2Q.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2Q::sTag()
   {
      return "Fwd::I2Q";
   }

   const std::size_t& I2Q::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2Q>(I2Q::sTag());
      return *i;
   }

   I2Q::I2Q()
      : IOperator(I2Q::sTag())
   {
   }

   I2Q::~I2Q()
   {
   }

}
}
}
