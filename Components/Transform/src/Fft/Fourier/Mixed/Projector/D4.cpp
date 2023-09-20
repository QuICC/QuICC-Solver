/**
 * @file D4.cpp
 * @brief Source of the implementation of the Fourier mixed D^4 projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D4.hpp"

// Project includes
//
#include "Types/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   D4::D4()
   {
   }

   D4::~D4()
   {
   }

   void D4::applyPreOperator(MatrixZ& out, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(out, in, 4, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
