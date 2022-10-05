/**
 * @file Laplh_1.cpp
 * @brief Source of the forward transform operator Forward::Laplh_1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh_1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh_1::sTag()
   {
      return "Fwd::Laplh_1";
   }

   std::string Laplh_1::sFormatted()
   {
      return "Forward::Laplh_1";
   }

   Laplh_1::Laplh_1()
      : IRegisterId<Laplh_1>(Laplh_1::sTag(), Laplh_1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
