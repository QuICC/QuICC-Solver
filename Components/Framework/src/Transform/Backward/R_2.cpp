/**
 * @file R_2.cpp
 * @brief Source of the backward transform operator Backard::R_2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_2::sTag()
   {
      return "Bwd::R_2";
   }

   std::string R_2::sFormatted()
   {
      return "Backard::R_2";
   }

   R_2::R_2()
      : IRegisterId<R_2>(R_2::sTag(), R_2::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
