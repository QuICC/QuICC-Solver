/**
 * @file D2.cpp
 * @brief Source of the implementation of the Fourier mixed D^2 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D2Base.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   void D2<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.outputDiff(rOut, 2, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
