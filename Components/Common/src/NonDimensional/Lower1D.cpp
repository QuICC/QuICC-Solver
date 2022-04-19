/**
 * @file Lower1D.cpp
 * @brief Source of the Lower1D nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Lower1D.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Lower1D::sTag()
   {
      return "lower1d";
   }

   std::string Lower1D::sFormatted()
   {
      return "Lower1D";
   }

   Lower1D::Lower1D(const MHDFloat value)
      : IRegisterId<Lower1D>(value, Lower1D::sTag(), Lower1D::sFormatted())
   {
   }

}
}
