/**
 * @file LaplhSin_1.cpp
 * @brief Source of the forward transform operator Forward::LaplhSin_1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/LaplhSin_1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string LaplhSin_1::sTag()
   {
      return "Fwd::LaplhSin_1";
   }

   std::string LaplhSin_1::sFormatted()
   {
      return "Forward::LaplhSin_1";
   }

   LaplhSin_1::LaplhSin_1()
      : IRegisterId<LaplhSin_1>(LaplhSin_1::sTag(), LaplhSin_1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
