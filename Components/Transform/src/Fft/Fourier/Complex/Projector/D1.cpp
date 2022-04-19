/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier complex D projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D1.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   D1::D1()
   {
   }

   D1::~D1()
   {
   }

   void D1::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(in, 1, this->mspSetup->boxScale());

      this->mBackend.output(rOut.data());
   }

   void D1::applyPostOperator(MatrixZ&) const
   {
   }

}
}
}
}
}
}
