/**
 * @file D1.cpp
 * @brief Source of the implementation of the Fourier complex D projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D1Base.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void D1<base_t>::applyPreOperator(MatrixZ& tmp, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->mBackend.inputDiff(tmp, in, 1, this->mspSetup->boxScale());
   }

}
}
}
}
}
}
