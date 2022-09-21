/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier mixed D integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   D1::D1()
   {
   }

   D1::~D1()
   {
   }

   void D1::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.outputDiff(rOut, 1, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
