/**
 * @file Entropy.cpp
 * @brief Source of the Entropy physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Entropy.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Entropy::sTag()
   {
      return "entropy";
   }

   std::string Entropy::sFormatted()
   {
      return "Entropy";
   }

   Entropy::Entropy()
      : IRegisterId<Entropy>(Entropy::sTag(), Entropy::sFormatted())
   {
   }

}
}
