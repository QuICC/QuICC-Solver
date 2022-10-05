/**
 * @file R_1.cpp
 * @brief Source of the backward transform operator Backard::R_1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1::sTag()
   {
      return "Bwd::R_1";
   }

   std::string R_1::sFormatted()
   {
      return "Backard::R_1";
   }

   R_1::R_1()
      : IRegisterId<R_1>(R_1::sTag(), R_1::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
