/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier complex D integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D1Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void D1<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.outputDiff(rOut, 1, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
