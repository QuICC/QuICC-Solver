/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier complex P projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/PBase.hpp"

// Pect includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void P<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mBackend.input(tmp, in);
   }

}
}
}
}
}
}
