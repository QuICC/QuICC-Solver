/** 
 * @file PowerD1.cpp
 * @brief Source of reductor projection operator PowerD1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerD1.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerD1::sTag()
   {
      return "Red::PowerD1";
   }

   const std::size_t& PowerD1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<PowerD1>(PowerD1::sTag());
      return *i;
   }

   PowerD1::PowerD1()
      : IOperator(PowerD1::sTag())
   {
   }

   PowerD1::~PowerD1()
   {
   }

}
}
}
