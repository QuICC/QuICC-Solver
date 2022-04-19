/**
 * @file Diagnostic.cpp
 * @brief Source of the Diagnostic PseudospectralTag
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PseudospectralTag/Diagnostic.hpp"

// Project includes
//

namespace QuICC {

namespace PseudospectralTag {

   std::string Diagnostic::sTag()
   {
      return "diagnostic";
   }

   std::string Diagnostic::sFormatted()
   {
      return "Diagnostic";
   }

   Diagnostic::Diagnostic()
      : IRegisterId<Diagnostic>(Diagnostic::sTag(), Diagnostic::sFormatted())
   {
   }

}
}
