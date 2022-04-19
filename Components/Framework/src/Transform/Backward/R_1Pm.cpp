/** 
 * @file R_1Pm.cpp
 * @brief Source of backward projection operator R_1Pm 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1Pm.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1Pm::sTag()
   {
      return "Bwd::R_1Pm";
   }

   const std::size_t& R_1Pm::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<R_1Pm>(R_1Pm::sTag());
      return *i;
   }

   R_1Pm::R_1Pm()
      : IOperator(R_1Pm::sTag())
   {
   }

   R_1Pm::~R_1Pm()
   {
   }

}
}
}
