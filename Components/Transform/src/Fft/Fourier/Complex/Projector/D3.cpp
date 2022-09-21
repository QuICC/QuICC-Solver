/**
 * @file D3.cpp
 * @brief Source of the implementation of the Fourier complex D^3 projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D3.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   D3::D3()
   {
   }

   D3::~D3()
   {
   }

   void D3::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(tmp, in, 3, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
