/**
 * @file P_Clean.cpp
 * @brief Source of the implementation of the Fourier complex P integrator, but 0 mode is cleaned
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/P_CleanBase.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void P_Clean<base_t>::initOperator() const
   {
      this->mBackend.initMeanBlocks(this->mspSetup->idBlocks());
   }

   void P_Clean<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.zeroMean(rOut);
      this->mBackend.output(rOut);
   }

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC
