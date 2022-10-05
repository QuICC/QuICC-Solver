/**
 * @file S.cpp
 * @brief Source of the forward transform operator Forward::S
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/S.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string S::sTag()
   {
      return "Fwd::S";
   }

   std::string S::sFormatted()
   {
      return "Forward::S";
   }

   S::S()
      : IRegisterId<S>(S::sTag(), S::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
