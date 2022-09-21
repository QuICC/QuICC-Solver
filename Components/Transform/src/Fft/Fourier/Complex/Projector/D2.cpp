/**
 * @file D2.cpp
 * @brief Source of the implementation of the Fourier complex D^2 projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   D2::D2()
   {
   }

   D2::~D2()
   {
   }

   void D2::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(tmp, in, 2, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
