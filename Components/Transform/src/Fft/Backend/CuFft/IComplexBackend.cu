/**
 * @file IComplexBackend.cu
 * @brief Source of the interface for a generic cuFFT based complex backend
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/IComplexBackend.hpp"

// Project includes
//
#include "Types/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   const int IComplexBackend::MIN_BATCH_BLOCKSIZE = 10;

   IComplexBackend::IComplexBackend()
      : IComplexBackend(1,1)
   {
   }

   IComplexBackend::IComplexBackend(const int nStreams, const int nBatches)
      : ICuFftBackend(nStreams, nBatches)
   {
   }

   IComplexBackend::~IComplexBackend()
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

   void IComplexBackend::init(const SetupType& setup) const
   {
      // Get size of positive and negative frequency parts
      int specSize = setup.specSize();
      this->mNegN = specSize/2;
      this->mPosN = this->mNegN + (specSize%2);
      this->mFwdSize = setup.fwdSize();
      this->mBwdSize = setup.bwdSize();

      // Recompute maximum number of batches
      this->mNBatches = std::min(setup.blockSize()/IComplexBackend::MIN_BATCH_BLOCKSIZE + 1, this->mNBatches);
      this->mNBatches = std::min(this->mNBatches, this->mNStreams);

      this->mBlockSize = setup.blockSize()/this->mNBatches;
   }

   void IComplexBackend::initMeanBlocks(const MatrixI& idBlocks) const
   {
      // Set the mean blocks
      int start = 0;
      for(int i = 0; i < idBlocks.rows(); ++i)
      {
         if(idBlocks(i,0) == 0)
         {
            this->mMeanBlocks.push_back(std::make_pair(start, idBlocks(i,1)));
         }

         start += idBlocks(i,1);
      }
   }

   void IComplexBackend::input(const MHDComplex* in) const
   {
      this->io(this->mTmp.data(), in);
   }

   void IComplexBackend::output(MHDComplex* out) const
   {
      this->io(out, this->mTmp.data());
   }

   void IComplexBackend::io(MHDComplex* out, const MHDComplex* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

   Array IComplexBackend::positiveK() const
   {
      return Array::LinSpaced(this->mPosN, 0, this->mPosN-1);
   }

   Array IComplexBackend::negativeK() const
   {
      return Array::LinSpaced(this->mNegN, 0, this->mNegN-1).array() - static_cast<double>(this->mNegN);
   }

}
}
}
}
}
