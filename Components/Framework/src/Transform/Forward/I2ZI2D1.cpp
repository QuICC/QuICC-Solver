/** 
 * @file I2ZI2D1.cpp
 * @brief Source of forward projection operator I2ZI2D1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2ZI2D1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2ZI2D1::sTag()
   {
      return "Fwd::I2ZI2D1";
   }

   const std::size_t& I2ZI2D1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2ZI2D1>(I2ZI2D1::sTag());
      return *i;
   }

   I2ZI2D1::I2ZI2D1()
      : IOperator(I2ZI2D1::sTag())
   {
   }

   I2ZI2D1::~I2ZI2D1()
   {
   }

}
}
}
