/**
 * @file ICuFftBackend.cu
 * @brief Source of the interface for a generic cuFFT based backend
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/ICuFftBackend.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Backend/CuFft/CheckCuda.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   ICuFftBackend::ICuFftBackend()
      : ICuFftBackend(1,1)
   {
   }

   ICuFftBackend::ICuFftBackend(const int nStreams, const int nBatches)
      : mNStreams(nStreams), mNBatches(nBatches)
   {
      assert(this->mNStreams > 0);

      this->initLibrary();
   }

   ICuFftBackend::~ICuFftBackend()
   {
      this->cleanupFft();
      this->cleanupLibrary();
   }

   void ICuFftBackend::initLibrary() const
   {
      // Init and register user with library
      Library::init();
      Library::registerUser();

      // Create CUDA stream(s)
      this->mStreams.reserve(this->mNStreams);
      for(int i = 0; i < this->mNStreams; i++)
      {
         cudaStream_t stream;
         CheckCuda(cudaStreamCreate(&stream), __LINE__);
         this->mStreams.push_back(stream);
      }
   }

   void ICuFftBackend::cleanupFft()
   {
       // Destroy plans
      if(this->mPlans.size() > 0)
      {
         for(auto it = this->mPlans.begin(); it != this->mPlans.end(); ++it)
         {
            cufftDestroy(*it);
         }

         for(auto it = this->mStreams.begin(); it != this->mStreams.end(); ++it)
         {
            CheckCuda(cudaStreamDestroy(*it), __LINE__);
         }
      }
   }

   void ICuFftBackend::synchronize() const
   {
      for(int j = 0; j < this->mNStreams; j++)
      {
         CheckCuda(cudaStreamSynchronize(this->mStreams.at(j)), __LINE__);
      }
   }

   void ICuFftBackend::synchronize(const int id) const
   {
      CheckCuda(cudaStreamSynchronize(this->mStreams.at(id)), __LINE__);
   }

   cudaStream_t ICuFftBackend::stream(const int id) const
   {
      return this->mStreams.at(id);
   }

   void ICuFftBackend::cleanupLibrary()
   {
      // Unregister user with library
      Library::unregisterUser();

      // cleanup library
      Library::cleanup();
   }

}
}
}
}
}
