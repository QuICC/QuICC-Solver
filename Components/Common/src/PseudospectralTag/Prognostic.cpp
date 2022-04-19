/**
 * @file Prognostic.cpp
 * @brief Source of the Prognostic PseudospectralTag
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PseudospectralTag/Prognostic.hpp"

// Project includes
//

namespace QuICC {

namespace PseudospectralTag {

   std::string Prognostic::sTag()
   {
      return "prognostic";
   }

   std::string Prognostic::sFormatted()
   {
      return "Prognostic";
   }

   Prognostic::Prognostic()
      : IRegisterId<Prognostic>(Prognostic::sTag(), Prognostic::sFormatted())
   {
   }

}
}
