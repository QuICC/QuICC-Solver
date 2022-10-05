/**
 * @file P.cpp
 * @brief Source of the forward transform operator Forward::P
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string P::sTag()
   {
      return "Fwd::P";
   }

   std::string P::sFormatted()
   {
      return "Forward::P";
   }

   P::P()
      : IRegisterId<P>(P::sTag(), P::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
