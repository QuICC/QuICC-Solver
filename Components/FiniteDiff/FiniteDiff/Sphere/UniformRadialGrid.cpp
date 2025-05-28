/**
 * @file UniformRadialGrid.cpp
 * @brief Setup uniform radial grid in sphere
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/UniformRadialGrid.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   void UniformRadialGrid::computeGrid(Internal::Array& igrid, const int size)
   {
      igrid.resize(size);
      Internal::MHDFloat delta = static_cast<Internal::MHDFloat>(1)/static_cast<Internal::MHDFloat>(size-1);
      for(int i = 0; i < size; i++)
      {
         igrid(i) = i*delta;
      }
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
