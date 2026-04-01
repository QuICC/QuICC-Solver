/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier mixed D projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D1Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   void D1<base_t>::applyPreOperator(MatrixZ& out, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->mBackend.inputDiff(out, in, 1, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
