/**
 * @file D4.cpp
 * @brief Source of the implementation of the Fourier mixed D^4 projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D4.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   void D4::applyPreOperator(MatrixZ& out, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->mBackend.inputDiff(out, in, 4, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
