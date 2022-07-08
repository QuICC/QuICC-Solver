/**
 * @file I2T.cpp
 * @brief Source of forward projection operator I2T
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2T.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2T::sTag()
   {
      return "Fwd::I2T";
   }

   const std::size_t& I2T::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2T>(I2T::sTag());
      return *i;
   }

   I2T::I2T()
      : IOperator(I2T::sTag())
   {
   }

   I2T::~I2T()
   {
   }

}
}
}
