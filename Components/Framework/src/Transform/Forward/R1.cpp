/** 
 * @file R1.cpp
 * @brief Source of forward projection operator R1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/R1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string R1::sTag()
   {
      return "Fwd::R1";
   }

   const std::size_t& R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<R1>(R1::sTag());
      return *i;
   }

   R1::R1()
      : IOperator(R1::sTag())
   {
   }

   R1::~R1()
   {
   }

}
}
}
