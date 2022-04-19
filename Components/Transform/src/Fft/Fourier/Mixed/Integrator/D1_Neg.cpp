/**
 * @file D1_Neg.cpp
 * @brief Source of the implementation of the Fourier mixed D integrator, but 0 mode is -P integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1_Neg.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   D1_Neg::D1_Neg()
   {
   }

   D1_Neg::~D1_Neg()
   {
   }

   void D1_Neg::applyPreOperator(MatrixZ& rOut, const Matrix& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void D1_Neg::applyPostOperator(MatrixZ& rOut) const
   {
      std::map<int,MHDComplex> mod = { {0, -1.0} };
      this->mBackend.outputDiff(rOut, 1, this->mspSetup->boxScale(), mod);
   }

}
}
}
}
}
}
