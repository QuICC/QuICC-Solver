/** 
 * @file S4.cpp
 * @brief Source of forward projection operator S4 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/S4.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string S4::sTag()
   {
      return "Fwd::S4";
   }

   const std::size_t& S4::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<S4>(S4::sTag());
      return *i;
   }

   S4::S4()
      : IOperator(S4::sTag())
   {
   }

   S4::~S4()
   {
   }

}
}
}
