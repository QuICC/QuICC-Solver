/**
 * @file Prandtl.cpp
 * @brief Source of the Prandtl nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Prandtl.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Prandtl::sTag()
   {
      return "prandtl";
   }

   std::string Prandtl::sFormatted()
   {
      return "Prandtl";
   }

   Prandtl::Prandtl(const MHDFloat value)
      : IRegisterId<Prandtl>(value, Prandtl::sTag(), Prandtl::sFormatted())
   {
   }

}
}
