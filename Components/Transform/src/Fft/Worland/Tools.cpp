/**
 * @file Tools.cpp
 * @brief Source of tools for Worland FFT transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "QuICC/Transform/Fft/Worland/Tools.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Tools {

   void computeGrid(Internal::Array& grid, const int size)
   {
      grid.resize(size);

      // Create Chebyshev grid for r = sqrt((x+1)/2)
      for(int k = 0; k < size; k++)
      {
         grid(k) = Internal::Math::cos(MHD_MP(0.5)*(Internal::Math::PI_long)*(Internal::MHDFloat(k)+ MHD_MP(0.5))/Internal::MHDFloat(size));
      }
   }

}
}
}
}
}
