/**
 * @file D2.cpp
 * @brief Source of the forward transform operator Forward::D2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/D2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string D2::sTag()
   {
      return "Fwd::D2";
   }

   std::string D2::sFormatted()
   {
      return "Forward::D2";
   }

   D2::D2()
      : IRegisterId<D2>(D2::sTag(), D2::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
