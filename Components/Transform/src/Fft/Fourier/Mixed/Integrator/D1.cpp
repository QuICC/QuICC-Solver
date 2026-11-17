/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier mixed D integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   void D1<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut) const
   {
      this->mBackend.outputDiff(rOut, 1, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
