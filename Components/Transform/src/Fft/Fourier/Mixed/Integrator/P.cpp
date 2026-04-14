/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier mixed P integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   void P<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut) const
   {
      this->mBackend.output(rOut);
   }

}
}
}
}
}
}
