/**
 * @file Mean.cpp
 * @brief Source of the implementation of the Fourier complex mean projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Mean.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   Mean::Mean()
   {
   }

   Mean::~Mean()
   {
   }

   void Mean::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.inputMean(tmp, in);
   }

}
}
}
}
}
}
