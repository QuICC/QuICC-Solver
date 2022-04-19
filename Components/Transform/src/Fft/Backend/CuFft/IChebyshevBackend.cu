/**
 * @file IChebyshevBackend.cu
 * @brief Source of the interface for a generic cuFFT based Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/IChebyshevBackend.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   const int IChebyshevBackend::MIN_BATCH_BLOCKSIZE = 10;

   IChebyshevBackend::IChebyshevBackend()
      : IChebyshevBackend(1,1)
   {
   }

   IChebyshevBackend::IChebyshevBackend(const int nStreams, const int nBatches)
      : ICuFftBackend(nStreams, nBatches)
   {
   }

   IChebyshevBackend::~IChebyshevBackend()
   {
      for(auto& d: this->mcuFwd)
      {
         cudaFree(d);
      }
      for(auto& d: this->mcuBwd)
      {
         cudaFree(d);
      }
   }

   void IChebyshevBackend::init(const SetupType& setup) const
   {
      this->mSpecSize = setup.specSize();
      this->mFwdSize = setup.fwdSize();
      this->mBwdSize = setup.bwdSize();

      if(this->mFwdSize%2 == 1)
      {
         throw std::logic_error("cuFFT backend only supports even resolutions");
      }

      // Recompute maximum number of batches
      this->mNBatches = std::min(setup.blockSize()/IChebyshevBackend::MIN_BATCH_BLOCKSIZE + 1, this->mNBatches);
      this->mNBatches = std::min(this->mNBatches, this->mNStreams);

      this->mBlockSize = setup.blockSize()/this->mNBatches;
   }

   void IChebyshevBackend::io(double* out, const double* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

}
}
}
}
}
