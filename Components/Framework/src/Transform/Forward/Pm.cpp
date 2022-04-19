/** 
 * @file Pm.cpp
 * @brief Source of forward projection operator Pm 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Pm.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Pm::sTag()
   {
      return "Fwd::Pm";
   }

   const std::size_t& Pm::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Pm>(Pm::sTag());
      return *i;
   }

   Pm::Pm()
      : IOperator(Pm::sTag())
   {
   }

   Pm::~Pm()
   {
   }

}
}
}
