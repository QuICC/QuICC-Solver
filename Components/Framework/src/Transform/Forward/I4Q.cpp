/**
 * @file I4Q.cpp
 * @brief Source of the forward transform operator Forward::I4Q
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4Q.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4Q::sTag()
   {
      return "Fwd::I4Q";
   }

   std::string I4Q::sFormatted()
   {
      return "Forward::I4Q";
   }

   I4Q::I4Q()
      : IRegisterId<I4Q>(I4Q::sTag(), I4Q::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
