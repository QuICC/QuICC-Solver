/**
 * @file Mean.cpp
 * @brief Source of the implementation of the Fourier complex mean integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/MeanBase.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void Mean<base_t>::initOperator() const
   {
      this->mBackend.initMeanBlocks(this->mspSetup->idBlocks());
   }

   void Mean<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut) const
   {
      this->mBackend.outputMean(rOut);
   }

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC
