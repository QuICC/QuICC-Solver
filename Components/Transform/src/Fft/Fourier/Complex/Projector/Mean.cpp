/**
 * @file Mean.cpp
 * @brief Source of the implementation of the Fourier complex mean projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/MeanBase.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void Mean<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputMean(tmp, in);
   }

}
}
}
}
}
}
