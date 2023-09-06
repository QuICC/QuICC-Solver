/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier mixed P projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/PBase.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   void P<base_t>::applyPreOperator(MatrixZ& out, const MatrixZ& in) const
   {
      this->mBackend.input(out, in);
   }
}
}
}
}
}
}
