/**
 * @file Trivial.cpp
 * @brief Source of the Trivial PseudospectralTag
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PseudospectralTag/Trivial.hpp"

// Project includes
//

namespace QuICC {

namespace PseudospectralTag {

   std::string Trivial::sTag()
   {
      return "trivial";
   }

   std::string Trivial::sFormatted()
   {
      return "Trivial";
   }

   Trivial::Trivial()
      : IRegisterId<Trivial>(Trivial::sTag(), Trivial::sFormatted())
   {
   }

}
}
