/** 
 * @file Q2.cpp
 * @brief Source of forward projection operator Q2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Q2.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Q2::sTag()
   {
      return "Fwd::Q2";
   }

   const std::size_t& Q2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Q2>(Q2::sTag());
      return *i;
   }

   Q2::Q2()
      : IOperator(Q2::sTag())
   {
   }

   Q2::~Q2()
   {
   }

}
}
}
