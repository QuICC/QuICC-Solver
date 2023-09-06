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
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D1Base.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   void D1<base_t>::applyPreOperator(MatrixZ& out, const MatrixZ& in) const
   {
      this->mBackend.inputDiff(out, in, 1, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
