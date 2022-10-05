/**
 * @file LaplhZR_1D1R1.cpp
 * @brief Source of the backward transform operator Backard::LaplhZR_1D1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/LaplhZR_1D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string LaplhZR_1D1R1::sTag()
   {
      return "Bwd::LaplhZR_1D1R1";
   }

   std::string LaplhZR_1D1R1::sFormatted()
   {
      return "Backard::LaplhZR_1D1R1";
   }

   LaplhZR_1D1R1::LaplhZR_1D1R1()
      : IRegisterId<LaplhZR_1D1R1>(LaplhZR_1D1R1::sTag(), LaplhZR_1D1R1::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
