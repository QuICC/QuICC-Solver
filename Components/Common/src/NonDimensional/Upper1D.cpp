/**
 * @file Upper1D.cpp
 * @brief Source of the Upper1D nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Upper1D.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Upper1D::sTag()
   {
      return "upper1d";
   }

   std::string Upper1D::sFormatted()
   {
      return "Upper1D";
   }

   Upper1D::Upper1D(const MHDFloat value)
      : IRegisterId<Upper1D>(value, Upper1D::sTag(), Upper1D::sFormatted())
   {
   }

}
}
