/** 
 * @file I4R_1Pm.cpp
 * @brief Source of forward projection operator I4R_1Pm 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4R_1Pm.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4R_1Pm::sTag()
   {
      return "Fwd::I4R_1Pm";
   }

   const std::size_t& I4R_1Pm::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4R_1Pm>(I4R_1Pm::sTag());
      return *i;
   }

   I4R_1Pm::I4R_1Pm()
      : IOperator(I4R_1Pm::sTag())
   {
   }

   I4R_1Pm::~I4R_1Pm()
   {
   }

}
}
}
