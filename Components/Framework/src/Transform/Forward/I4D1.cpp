/**
 * @file I4D1.cpp
 * @brief Source of the forward transform operator Forward::I4D1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4D1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4D1::sTag()
   {
      return "Fwd::I4D1";
   }

   std::string I4D1::sFormatted()
   {
      return "Forward::I4D1";
   }

   I4D1::I4D1()
      : IRegisterId<I4D1>(I4D1::sTag(), I4D1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
