/** 
 * @file PowerR2.cpp
 * @brief Source of reductor projection operator PowerR2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerR2.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerR2::sTag()
   {
      return "Red::PowerR2";
   }

   const std::size_t& PowerR2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<PowerR2>(PowerR2::sTag());
      return *i;
   }

   PowerR2::PowerR2()
      : IOperator(PowerR2::sTag())
   {
   }

   PowerR2::~PowerR2()
   {
   }

}
}
}
