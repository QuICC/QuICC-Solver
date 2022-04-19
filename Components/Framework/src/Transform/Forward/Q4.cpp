/** 
 * @file Q4.cpp
 * @brief Source of forward projection operator Q4 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Q4.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Q4::sTag()
   {
      return "Fwd::Q4";
   }

   const std::size_t& Q4::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Q4>(Q4::sTag());
      return *i;
   }

   Q4::Q4()
      : IOperator(Q4::sTag())
   {
   }

   Q4::~Q4()
   {
   }

}
}
}
