/**
 * @file I4P.cpp
 * @brief Source of the forward transform operator Forward::I4P
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4P::sTag()
   {
      return "Fwd::I4P";
   }

   std::string I4P::sFormatted()
   {
      return "Forward::I4P";
   }

   I4P::I4P()
      : IRegisterId<I4P>(I4P::sTag(), I4P::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
