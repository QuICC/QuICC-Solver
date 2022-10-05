/**
 * @file R_1Pm.cpp
 * @brief Source of the backward transform operator Backard::R_1Pm
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1Pm.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1Pm::sTag()
   {
      return "Bwd::R_1Pm";
   }

   std::string R_1Pm::sFormatted()
   {
      return "Backard::R_1Pm";
   }

   R_1Pm::R_1Pm()
      : IRegisterId<R_1Pm>(R_1Pm::sTag(), R_1Pm::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
