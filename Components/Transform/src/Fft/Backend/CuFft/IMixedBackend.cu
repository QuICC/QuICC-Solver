/**
 * @file IMixedBackend.cu
 * @brief Source of the interface for a generic cuFFT based mixed backend
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/IMixedBackend.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   const int IMixedBackend::MIN_BATCH_BLOCKSIZE = 10;

   IMixedBackend::IMixedBackend()
      : IMixedBackend(1,1)
   {
   }

   IMixedBackend::IMixedBackend(const int nStreams, const int nBatches)
      : ICuFftBackend(nStreams, nBatches)
   {
   }

   IMixedBackend::~IMixedBackend()
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

   void IMixedBackend::init(const SetupType& setup) const
   {
      this->mSpecSize = setup.specSize();
      this->mFwdSize = setup.fwdSize();
      this->mBwdSize = setup.bwdSize();

      // Recompute maximum number of batches
      this->mNBatches = std::min(setup.blockSize()/IMixedBackend::MIN_BATCH_BLOCKSIZE + 1, this->mNBatches);
      this->mNBatches = std::min(this->mNBatches, this->mNStreams);

      this->mBlockSize = setup.blockSize()/this->mNBatches;
   }

   Array IMixedBackend::positiveK() const
   {
      return Array::LinSpaced(this->mSpecSize, 0, this->mSpecSize-1);
   }

}
}
}
}
}
