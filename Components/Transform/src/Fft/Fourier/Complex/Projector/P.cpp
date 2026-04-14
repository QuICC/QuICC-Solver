/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier complex P projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/PBase.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void P<base_t>::applyPreOperator(MatrixZ& tmp, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->mBackend.input(tmp, in);
   }

}
}
}
}
}
}
