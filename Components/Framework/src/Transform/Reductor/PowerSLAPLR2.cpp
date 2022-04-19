/** 
 * @file PowerSLAPLR2.cpp
 * @brief Source of reductor projection operator PowerSLAPLR2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerSLAPLR2.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerSLAPLR2::sTag()
   {
      return "Red::PowerSLAPLR2";
   }

   const std::size_t& PowerSLAPLR2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<PowerSLAPLR2>(PowerSLAPLR2::sTag());
      return *i;
   }

   PowerSLAPLR2::PowerSLAPLR2()
      : IOperator(PowerSLAPLR2::sTag())
   {
   }

   PowerSLAPLR2::~PowerSLAPLR2()
   {
   }

}
}
}
