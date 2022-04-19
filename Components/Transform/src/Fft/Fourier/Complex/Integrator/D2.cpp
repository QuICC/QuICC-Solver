/**
 * @file D2.cpp
 * @brief Source of the implementation of the Fourier complex D^2 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   D2::D2()
   {
   }

   D2::~D2()
   {
   }

   void D2::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void D2::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.outputDiff(rOut, 2, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
