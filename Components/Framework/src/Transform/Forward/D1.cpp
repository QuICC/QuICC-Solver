/**
 * @file D1.cpp
 * @brief Source of the forward transform operator Forward::D1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/D1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string D1::sTag()
   {
      return "Fwd::D1";
   }

   std::string D1::sFormatted()
   {
      return "Forward::D1";
   }

   D1::D1()
      : IRegisterId<D1>(D1::sTag(), D1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
