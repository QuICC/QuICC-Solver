/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier mixed P projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   P::P()
   {
   }

   P::~P()
   {
   }

   void P::applyPreOperator(MatrixZ& out, const MatrixZ& in) const
   {
      this->mBackend.input(out, in);
   }

}
}
}
}
}
}
