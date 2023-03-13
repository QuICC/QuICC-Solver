/**
 * @file D3.cpp
 * @brief Source of the implementation of the Fourier complex D^3 projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D3Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void D3<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(tmp, in, 3, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
