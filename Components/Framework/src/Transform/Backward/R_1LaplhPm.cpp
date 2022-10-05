/**
 * @file R_1LaplhPm.cpp
 * @brief Source of the backward transform operator Backard::R_1LaplhPm
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1LaplhPm.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1LaplhPm::sTag()
   {
      return "Bwd::R_1LaplhPm";
   }

   std::string R_1LaplhPm::sFormatted()
   {
      return "Backard::R_1LaplhPm";
   }

   R_1LaplhPm::R_1LaplhPm()
      : IRegisterId<R_1LaplhPm>(R_1LaplhPm::sTag(), R_1LaplhPm::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
