/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier complex P integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/PBase.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void P<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.output(rOut);
   }

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC
