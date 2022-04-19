/** 
 * @file S2.cpp
 * @brief Source of forward projection operator S2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/S2.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string S2::sTag()
   {
      return "Fwd::S2";
   }

   const std::size_t& S2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<S2>(S2::sTag());
      return *i;
   }

   S2::S2()
      : IOperator(S2::sTag())
   {
   }

   S2::~S2()
   {
   }

}
}
}
