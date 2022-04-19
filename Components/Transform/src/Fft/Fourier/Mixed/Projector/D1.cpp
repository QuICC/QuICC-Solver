/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier mixed D projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   D1::D1()
   {
   }

   D1::~D1()
   {
   }

   void D1::applyPreOperator(Matrix& rOut, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(in, 1, this->mspSetup->boxScale());

      this->mBackend.output(rOut.data());
   }

   void D1::applyPostOperator(Matrix&) const
   {
   }

}
}
}
}
}
}
