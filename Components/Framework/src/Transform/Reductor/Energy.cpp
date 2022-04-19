/** 
 * @file Energy.cpp
 * @brief Source of reductor projection operator Energy 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/Energy.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string Energy::sTag()
   {
      return "Red::Energy";
   }

   const std::size_t& Energy::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Energy>(Energy::sTag());
      return *i;
   }

   Energy::Energy()
      : IOperator(Energy::sTag())
   {
   }

   Energy::~Energy()
   {
   }

}
}
}
