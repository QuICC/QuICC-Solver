/** 
 * @file Pol.cpp
 * @brief Source of forward projection operator Pol 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Pol.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Pol::sTag()
   {
      return "Fwd::Pol";
   }

   const std::size_t& Pol::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Pol>(Pol::sTag());
      return *i;
   }

   Pol::Pol()
      : IOperator(Pol::sTag())
   {
   }

   Pol::~Pol()
   {
   }

}
}
}
