/**
 * @file Pm.cpp
 * @brief Source of the forward transform operator Forward::Pm
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Pm.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Pm::sTag()
   {
      return "Fwd::Pm";
   }

   std::string Pm::sFormatted()
   {
      return "Forward::Pm";
   }

   Pm::Pm()
      : IRegisterId<Pm>(Pm::sTag(), Pm::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
