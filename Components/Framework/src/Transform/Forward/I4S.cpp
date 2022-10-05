/**
 * @file I4S.cpp
 * @brief Source of the forward transform operator Forward::I4S
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4S.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4S::sTag()
   {
      return "Fwd::I4S";
   }

   std::string I4S::sFormatted()
   {
      return "Forward::I4S";
   }

   I4S::I4S()
      : IRegisterId<I4S>(I4S::sTag(), I4S::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
