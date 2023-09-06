/**
 * @file D1_Neg.cpp
 * @brief Source of the implementation of the Fourier complex D integrator, but 0 mode is -P integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D1_NegBase.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void D1_Neg<base_t>::initOperator() const
   {
      this->mBackend.initMeanBlocks(this->mspSetup->idBlocks());
   }

   void D1_Neg<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.extractMean(rOut);

      this->mBackend.outputDiff(rOut, 1, this->mspSetup->boxScale());

      this->mBackend.setMean(rOut, -1.0);
   }

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC
