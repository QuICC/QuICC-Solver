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
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/P.hpp"

// Pect includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   P::P()
   {
   }

   P::~P()
   {
   }

   void P::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      this->mBackend.input(in);

      this->mBackend.output(rOut.data());
   }

   void P::applyPostOperator(MatrixZ&) const
   {
   }

}
}
}
}
}
}
