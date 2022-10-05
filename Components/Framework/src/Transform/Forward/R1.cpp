/**
 * @file R1.cpp
 * @brief Source of the forward transform operator Forward::R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string R1::sTag()
   {
      return "Fwd::R1";
   }

   std::string R1::sFormatted()
   {
      return "Forward::R1";
   }

   R1::R1()
      : IRegisterId<R1>(R1::sTag(), R1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
