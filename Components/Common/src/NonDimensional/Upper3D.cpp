/**
 * @file Upper3D.cpp
 * @brief Source of the Upper3D nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Upper3D.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Upper3D::sTag()
   {
      return "upper3d";
   }

   std::string Upper3D::sFormatted()
   {
      return "Upper3D";
   }

   Upper3D::Upper3D(const MHDFloat value)
      : IRegisterId<Upper3D>(value, Upper3D::sTag(), Upper3D::sFormatted())
   {
   }

}
}
