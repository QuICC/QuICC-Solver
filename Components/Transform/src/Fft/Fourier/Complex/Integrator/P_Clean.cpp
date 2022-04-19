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
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/P_Clean.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   P_Clean::P_Clean()
   {
   }

   P_Clean::~P_Clean()
   {
   }

   void P_Clean::initOperator() const
   {
      this->mBackend.initMeanBlocks(this->mspSetup->idBlocks());
   }

   void P_Clean::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void P_Clean::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.outputZeroMean(rOut);
   }

}
}
}
}
}
}
