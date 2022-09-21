/**
 * @file D3.cpp
 * @brief Source of the implementation of the Fourier mixed D^3 projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D3.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   D3::D3()
   {
   }

   D3::~D3()
   {
   }

   void D3::applyPreOperator(MatrixZ& out, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(out, in, 3, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
