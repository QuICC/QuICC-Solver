/**
 * @file R_1D1R1.cpp
 * @brief Source of the backward transform operator Backard::R_1D1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1D1R1::sTag()
   {
      return "Bwd::R_1D1R1";
   }

   std::string R_1D1R1::sFormatted()
   {
      return "Backard::R_1D1R1";
   }

   R_1D1R1::R_1D1R1()
      : IRegisterId<R_1D1R1>(R_1D1R1::sTag(), R_1D1R1::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
