/**
 * @file D1_Neg.cpp
 * @brief Source of the implementation of the Fourier complex D integrator, but 0 mode is -P integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D1_Neg.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   D1_Neg::D1_Neg()
   {
   }

   D1_Neg::~D1_Neg()
   {
   }

   void D1_Neg::initOperator() const
   {
      this->mBackend.initMeanBlocks(this->mspSetup->idBlocks());
   }

   void D1_Neg::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void D1_Neg::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.extractMean();

      this->mBackend.outputDiff(rOut, 1, this->mspSetup->boxScale());

      this->mBackend.setMean(rOut, -1.0);
   }

}
}
}
}
}
}
