/**
 * @file Sin_1LaplhDphi.cpp
 * @brief Source of the backward transform operator Backard::Sin_1LaplhDphi
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1LaplhDphi.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1LaplhDphi::sTag()
   {
      return "Bwd::Sin_1LaplhDphi";
   }

   std::string Sin_1LaplhDphi::sFormatted()
   {
      return "Backard::Sin_1LaplhDphi";
   }

   Sin_1LaplhDphi::Sin_1LaplhDphi()
      : IRegisterId<Sin_1LaplhDphi>(Sin_1LaplhDphi::sTag(), Sin_1LaplhDphi::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
