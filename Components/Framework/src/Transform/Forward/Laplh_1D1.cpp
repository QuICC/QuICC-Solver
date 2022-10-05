/**
 * @file Laplh_1D1.cpp
 * @brief Source of the forward transform operator Forward::Laplh_1D1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh_1D1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh_1D1::sTag()
   {
      return "Fwd::Laplh_1D1";
   }

   std::string Laplh_1D1::sFormatted()
   {
      return "Forward::Laplh_1D1";
   }

   Laplh_1D1::Laplh_1D1()
      : IRegisterId<Laplh_1D1>(Laplh_1D1::sTag(), Laplh_1D1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
