/** 
 * @file RadialPowerDivR1D1R1.cpp
 * @brief Source of reductor projection operator RadialPowerDivR1D1R1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/RadialPowerDivR1D1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string RadialPowerDivR1D1R1::sTag()
   {
      return "Red::RadialPowerDivR1D1R1";
   }

   const std::size_t& RadialPowerDivR1D1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<RadialPowerDivR1D1R1>(RadialPowerDivR1D1R1::sTag());
      return *i;
   }

   RadialPowerDivR1D1R1::RadialPowerDivR1D1R1()
      : IOperator(RadialPowerDivR1D1R1::sTag())
   {
   }

   RadialPowerDivR1D1R1::~RadialPowerDivR1D1R1()
   {
   }

}
}
}
