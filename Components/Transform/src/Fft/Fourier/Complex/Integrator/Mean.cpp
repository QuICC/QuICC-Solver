/**
 * @file Mean.cpp
 * @brief Source of the implementation of the Fourier complex mean integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/Mean.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   Mean::Mean()
   {
   }

   Mean::~Mean()
   {
   }

   void Mean::initOperator() const
   {
      this->mBackend.initMeanBlocks(this->mspSetup->idBlocks());
   }

   void Mean::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void Mean::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.outputMean(rOut);
   }

}
}
}
}
}
}
