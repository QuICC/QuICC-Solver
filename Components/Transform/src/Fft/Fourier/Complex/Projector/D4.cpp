/**
 * @file D4.cpp
 * @brief Source of the implementation of the Fourier complex D^4 projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D4Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void D4<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(tmp, in, 4, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
