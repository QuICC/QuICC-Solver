/**
 * @file D2.cpp
 * @brief Source of the implementation of the Fourier mixed D^2 projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   D2::D2()
   {
   }

   D2::~D2()
   {
   }

   void D2::applyPreOperator(MatrixZ& out, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(out, in, 2, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
