/**
 * @file Lower3D.cpp
 * @brief Source of the Lower3D nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Lower3D.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Lower3D::sTag()
   {
      return "lower3d";
   }

   std::string Lower3D::sFormatted()
   {
      return "Lower3D";
   }

   Lower3D::Lower3D(const MHDFloat value)
      : IRegisterId<Lower3D>(value, Lower3D::sTag(), Lower3D::sFormatted())
   {
   }

}
}
