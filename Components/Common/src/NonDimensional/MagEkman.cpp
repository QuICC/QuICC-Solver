/**
 * @file MagEkman.cpp
 * @brief Source of the magnetic Ekman nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/MagEkman.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string MagEkman::sTag()
   {
      return "magnetic_ekman";
   }

   std::string MagEkman::sFormatted()
   {
      return "magnetic Ekman";
   }

   MagEkman::MagEkman(const MHDFloat value)
      : IRegisterId<MagEkman>(value, MagEkman::sTag(), MagEkman::sFormatted())
   {
   }

}
}
