/** 
 * @file RadialPower.cpp
 * @brief Source of reductor projection operator RadialPower 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/RadialPower.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string RadialPower::sTag()
   {
      return "Red::RadialPower";
   }

   const std::size_t& RadialPower::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<RadialPower>(RadialPower::sTag());
      return *i;
   }

   RadialPower::RadialPower()
      : IOperator(RadialPower::sTag())
   {
   }

   RadialPower::~RadialPower()
   {
   }

}
}
}
