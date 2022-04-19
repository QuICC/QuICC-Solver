/**
 * @file Stencil.cpp
 * @brief Source of the Stencil ModelOperatorBoundary
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperatorBoundary {

   std::string Stencil::sTag()
   {
      return "stencil";
   }

   std::string Stencil::sFormatted()
   {
      return "Stencil";
   }

   Stencil::Stencil()
      : IRegisterId<Stencil>(Stencil::sTag(), Stencil::sFormatted())
   {
   }

}
}
