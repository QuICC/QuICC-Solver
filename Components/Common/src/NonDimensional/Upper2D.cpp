/**
 * @file Upper2D.cpp
 * @brief Source of the Upper2D nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Upper2D.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Upper2D::sTag()
   {
      return "upper2d";
   }

   std::string Upper2D::sFormatted()
   {
      return "Upper2D";
   }

   Upper2D::Upper2D(const MHDFloat value)
      : IRegisterId<Upper2D>(value, Upper2D::sTag(), Upper2D::sFormatted())
   {
   }

}
}
