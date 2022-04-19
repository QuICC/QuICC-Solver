/**
 * @file Wrapper.cpp
 * @brief Source of the Wrapper PseudospectralTag
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PseudospectralTag/Wrapper.hpp"

// Project includes
//

namespace QuICC {

namespace PseudospectralTag {

   std::string Wrapper::sTag()
   {
      return "wrapper";
   }

   std::string Wrapper::sFormatted()
   {
      return "Wrapper";
   }

   Wrapper::Wrapper()
      : IRegisterId<Wrapper>(Wrapper::sTag(), Wrapper::sFormatted())
   {
   }

}
}
