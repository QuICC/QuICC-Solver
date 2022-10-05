/**
 * @file T.cpp
 * @brief Source of the forward transform operator Forward::T
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/T.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string T::sTag()
   {
      return "Fwd::T";
   }

   std::string T::sFormatted()
   {
      return "Forward::T";
   }

   T::T()
      : IRegisterId<T>(T::sTag(), T::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
