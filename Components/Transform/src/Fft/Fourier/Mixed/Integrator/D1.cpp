/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier mixed D integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1Base.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

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
