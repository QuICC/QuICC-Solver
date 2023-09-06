/**
 * @file D2.cpp
 * @brief Source of the implementation of the Fourier complex D^2 integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D2Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void D2<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.outputDiff(rOut, 2, this->mspSetup->boxScale());
   }

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC
