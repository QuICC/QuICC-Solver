/**
 * @file I6R_1D1R1ZI4.cpp
 * @brief Source of the forward transform operator Forward::I6R_1D1R1ZI4
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6R_1D1R1ZI4.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6R_1D1R1ZI4::sTag()
   {
      return "Fwd::I6R_1D1R1ZI4";
   }

   std::string I6R_1D1R1ZI4::sFormatted()
   {
      return "Forward::I6R_1D1R1ZI4";
   }

   I6R_1D1R1ZI4::I6R_1D1R1ZI4()
      : IRegisterId<I6R_1D1R1ZI4>(I6R_1D1R1ZI4::sTag(), I6R_1D1R1ZI4::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
