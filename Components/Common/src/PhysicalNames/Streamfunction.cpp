/**
 * @file Streamfunction.cpp
 * @brief Source of the Streamfunction physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Streamfunction.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Streamfunction::sTag()
   {
      return "streamfunction";
   }

   std::string Streamfunction::sFormatted()
   {
      return "Streamfunction";
   }

   Streamfunction::Streamfunction()
      : IRegisterId<Streamfunction>(Streamfunction::sTag(), Streamfunction::sFormatted())
   {
   }

}
}
