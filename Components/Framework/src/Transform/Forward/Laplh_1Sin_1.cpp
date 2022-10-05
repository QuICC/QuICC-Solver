/**
 * @file Laplh_1Sin_1.cpp
 * @brief Source of the forward transform operator Forward::Laplh_1Sin_1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh_1Sin_1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh_1Sin_1::sTag()
   {
      return "Fwd::Laplh_1Sin_1";
   }

   std::string Laplh_1Sin_1::sFormatted()
   {
      return "Forward::Laplh_1Sin_1";
   }

   Laplh_1Sin_1::Laplh_1Sin_1()
      : IRegisterId<Laplh_1Sin_1>(Laplh_1Sin_1::sTag(), Laplh_1Sin_1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
