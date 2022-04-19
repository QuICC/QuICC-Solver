/** 
 * @file EnergySLAPLR2.cpp
 * @brief Source of reductor projection operator EnergySLAPLR2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergySLAPLR2.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergySLAPLR2::sTag()
   {
      return "Red::EnergySLAPLR2";
   }

   const std::size_t& EnergySLAPLR2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<EnergySLAPLR2>(EnergySLAPLR2::sTag());
      return *i;
   }

   EnergySLAPLR2::EnergySLAPLR2()
      : IOperator(EnergySLAPLR2::sTag())
   {
   }

   EnergySLAPLR2::~EnergySLAPLR2()
   {
   }

}
}
}
