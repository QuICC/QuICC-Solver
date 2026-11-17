/**
 * @file D1_P.cpp
 * @brief Source of the implementation of the Fourier mixed D integrator, but 0 mode is P integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1_P.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   void D1_P<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut) const
   {
      std::map<int,MHDComplex> mod = { {0, 1.0} };
      this->mBackend.outputDiff(rOut, 1, this->mspSetup->boxScale(), mod);
   }

}
}
}
}
}
}
