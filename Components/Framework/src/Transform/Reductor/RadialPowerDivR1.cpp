/** 
 * @file RadialPowerDivR1.cpp
 * @brief Source of reductor projection operator RadialPowerDivR1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/RadialPowerDivR1.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string RadialPowerDivR1::sTag()
   {
      return "Red::RadialPowerDivR1";
   }

   const std::size_t& RadialPowerDivR1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<RadialPowerDivR1>(RadialPowerDivR1::sTag());
      return *i;
   }

   RadialPowerDivR1::RadialPowerDivR1()
      : IOperator(RadialPowerDivR1::sTag())
   {
   }

   RadialPowerDivR1::~RadialPowerDivR1()
   {
   }

}
}
}
