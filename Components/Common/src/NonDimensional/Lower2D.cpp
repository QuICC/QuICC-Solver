/**
 * @file Lower2D.cpp
 * @brief Source of the Lower2D nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Lower2D.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Lower2D::sTag()
   {
      return "lower2d";
   }

   std::string Lower2D::sFormatted()
   {
      return "Lower2D";
   }

   Lower2D::Lower2D(const MHDFloat value)
      : IRegisterId<Lower2D>(value, Lower2D::sTag(), Lower2D::sFormatted())
   {
   }

}
}
