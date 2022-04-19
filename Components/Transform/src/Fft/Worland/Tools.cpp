/**
 * @file Tools.cpp
 * @brief Source of tools for Worland FFT transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Tools.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Tools {

   void computeGrid(internal::Array& grid, const int size)
   {
      grid.resize(size);

      // Create Chebyshev grid for r = sqrt((x+1)/2)
      for(int k = 0; k < size; k++)
      {
         grid(k) = precision::cos(MHD_MP(0.5)*(Precision::PI_long)*(internal::MHDFloat(k)+ MHD_MP(0.5))/internal::MHDFloat(size));
      }
   }

}
}
}
}
}
