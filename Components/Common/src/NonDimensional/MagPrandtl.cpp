/**
 * @file MagPrandtl.cpp
 * @brief Source of the magnetic Prandtl nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/MagPrandtl.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string MagPrandtl::sTag()
   {
      return "magnetic_prandtl";
   }

   std::string MagPrandtl::sFormatted()
   {
      return "magnetic Prandtl";
   }

   MagPrandtl::MagPrandtl(const MHDFloat value)
      : IRegisterId<MagPrandtl>(value, MagPrandtl::sTag(), MagPrandtl::sFormatted())
   {
   }

}
}
