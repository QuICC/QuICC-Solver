/**
 * @file IChebyshevBackend.cpp
 * @brief Source of the interface for a generic FFTW based Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/IChebyshevBackend.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   IChebyshevBackend::IChebyshevBackend()
   {
   }

   IChebyshevBackend::~IChebyshevBackend()
   {
   }

   void IChebyshevBackend::init(const SetupType& setup) const
   {
      this->mSpecSize = setup.specSize();
      this->mBlockSize = setup.blockSize();
   }

   void IChebyshevBackend::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

}
}
}
}
}
