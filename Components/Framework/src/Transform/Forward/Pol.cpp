/**
 * @file Pol.cpp
 * @brief Source of the forward transform operator Forward::Pol
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Pol.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Pol::sTag()
   {
      return "Fwd::Pol";
   }

   std::string Pol::sFormatted()
   {
      return "Forward::Pol";
   }

   Pol::Pol()
      : IRegisterId<Pol>(Pol::sTag(), Pol::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
