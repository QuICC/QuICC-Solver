/**
 * @file MagReynolds.cpp
 * @brief Source of the magnetic Reynolds nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/MagReynolds.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string MagReynolds::sTag()
   {
      return "magnetic_reynolds";
   }

   std::string MagReynolds::sFormatted()
   {
      return "magnetic Reynolds";
   }

   MagReynolds::MagReynolds(const MHDFloat value)
      : IRegisterId<MagReynolds>(value, MagReynolds::sTag(), MagReynolds::sFormatted())
   {
   }

}
}
