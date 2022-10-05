/**
 * @file I6Laplh.cpp
 * @brief Source of the forward transform operator Forward::I6Laplh
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6Laplh.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6Laplh::sTag()
   {
      return "Fwd::I6Laplh";
   }

   std::string I6Laplh::sFormatted()
   {
      return "Forward::I6Laplh";
   }

   I6Laplh::I6Laplh()
      : IRegisterId<I6Laplh>(I6Laplh::sTag(), I6Laplh::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
