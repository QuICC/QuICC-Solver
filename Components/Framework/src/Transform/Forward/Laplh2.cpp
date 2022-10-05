/**
 * @file Laplh2.cpp
 * @brief Source of the forward transform operator Forward::Laplh2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh2::sTag()
   {
      return "Fwd::Laplh2";
   }

   std::string Laplh2::sFormatted()
   {
      return "Forward::Laplh2";
   }

   Laplh2::Laplh2()
      : IRegisterId<Laplh2>(Laplh2::sTag(), Laplh2::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
