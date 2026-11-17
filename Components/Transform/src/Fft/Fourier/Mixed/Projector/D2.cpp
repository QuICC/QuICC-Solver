/**
 * @file D2.cpp
 * @brief Source of the implementation of the Fourier mixed D^2 projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D2Base.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   void D2<base_t>::applyPreOperator(MatrixZ& out, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->mBackend.inputDiff(out, in, 2, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
