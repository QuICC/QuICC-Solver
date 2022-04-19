/** 
 * @file I2D1.cpp
 * @brief Source of forward projection operator I2D1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2D1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2D1::sTag()
   {
      return "Fwd::I2D1";
   }

   const std::size_t& I2D1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2D1>(I2D1::sTag());
      return *i;
   }

   I2D1::I2D1()
      : IOperator(I2D1::sTag())
   {
   }

   I2D1::~I2D1()
   {
   }

}
}
}
