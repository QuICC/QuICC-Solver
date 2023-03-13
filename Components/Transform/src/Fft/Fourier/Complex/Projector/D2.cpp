/**
 * @file D2.cpp
 * @brief Source of the implementation of the Fourier complex D^2 projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D2Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void D2<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(tmp, in, 2, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
