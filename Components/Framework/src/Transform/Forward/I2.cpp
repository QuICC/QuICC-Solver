/** 
 * @file I2.cpp
 * @brief Source of forward projection operator I2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2::sTag()
   {
      return "Fwd::I2";
   }

   const std::size_t& I2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2>(I2::sTag());
      return *i;
   }

   I2::I2()
      : IOperator(I2::sTag())
   {
   }

   I2::~I2()
   {
   }

}
}
}
