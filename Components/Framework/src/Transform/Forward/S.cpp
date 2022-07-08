/** 
 * @file S.cpp
 * @brief Source of forward projection operator S 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/S.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string S::sTag()
   {
      return "Fwd::S";
   }

   const std::size_t& S::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<S>(S::sTag());
      return *i;
   }

   S::S()
      : IOperator(S::sTag())
   {
   }

   S::~S()
   {
   }

}
}
}
