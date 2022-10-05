/**
 * @file Laplh.cpp
 * @brief Source of the forward transform operator Forward::Laplh
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh::sTag()
   {
      return "Fwd::Laplh";
   }

   std::string Laplh::sFormatted()
   {
      return "Forward::Laplh";
   }

   Laplh::Laplh()
      : IRegisterId<Laplh>(Laplh::sTag(), Laplh::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
