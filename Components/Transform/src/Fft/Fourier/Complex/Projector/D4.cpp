/**
 * @file D4.cpp
 * @brief Source of the implementation of the Fourier complex D^4 projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D4.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   D4::D4()
   {
   }

   D4::~D4()
   {
   }

   void D4::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(tmp, in, 4, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
