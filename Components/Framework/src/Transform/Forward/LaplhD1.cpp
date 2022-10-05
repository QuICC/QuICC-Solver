/**
 * @file LaplhD1.cpp
 * @brief Source of the forward transform operator Forward::LaplhD1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/LaplhD1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string LaplhD1::sTag()
   {
      return "Fwd::LaplhD1";
   }

   std::string LaplhD1::sFormatted()
   {
      return "Forward::LaplhD1";
   }

   LaplhD1::LaplhD1()
      : IRegisterId<LaplhD1>(LaplhD1::sTag(), LaplhD1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
