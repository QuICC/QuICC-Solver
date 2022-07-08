/** 
 * @file Q.cpp
 * @brief Source of forward projection operator Q 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Q.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Q::sTag()
   {
      return "Fwd::Q";
   }

   const std::size_t& Q::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Q>(Q::sTag());
      return *i;
   }

   Q::Q()
      : IOperator(Q::sTag())
   {
   }

   Q::~Q()
   {
   }

}
}
}
