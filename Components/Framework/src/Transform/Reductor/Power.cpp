/** 
 * @file Power.cpp
 * @brief Source of reductor projection operator Power 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/Power.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string Power::sTag()
   {
      return "Red::Power";
   }

   const std::size_t& Power::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Power>(Power::sTag());
      return *i;
   }

   Power::Power()
      : IOperator(Power::sTag())
   {
   }

   Power::~Power()
   {
   }

}
}
}
