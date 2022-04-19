/** 
 * @file PowerD1R1.cpp
 * @brief Source of reductor projection operator PowerD1R1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerD1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerD1R1::sTag()
   {
      return "Red::PowerD1R1";
   }

   const std::size_t& PowerD1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<PowerD1R1>(PowerD1R1::sTag());
      return *i;
   }

   PowerD1R1::PowerD1R1()
      : IOperator(PowerD1R1::sTag())
   {
   }

   PowerD1R1::~PowerD1R1()
   {
   }

}
}
}
