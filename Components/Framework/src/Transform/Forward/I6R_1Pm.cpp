/** 
 * @file I6R_1Pm.cpp
 * @brief Source of forward projection operator I6R_1Pm 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6R_1Pm.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6R_1Pm::sTag()
   {
      return "Fwd::I6R_1Pm";
   }

   const std::size_t& I6R_1Pm::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I6R_1Pm>(I6R_1Pm::sTag());
      return *i;
   }

   I6R_1Pm::I6R_1Pm()
      : IOperator(I6R_1Pm::sTag())
   {
   }

   I6R_1Pm::~I6R_1Pm()
   {
   }

}
}
}
