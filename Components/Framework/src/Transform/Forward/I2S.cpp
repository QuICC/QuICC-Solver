/**
 * @file I2S.cpp
 * @brief Source of forward projection operator I2S
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2S.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2S::sTag()
   {
      return "Fwd::I2S";
   }

   const std::size_t& I2S::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2S>(I2S::sTag());
      return *i;
   }

   I2S::I2S()
      : IOperator(I2S::sTag())
   {
   }

   I2S::~I2S()
   {
   }

}
}
}
