/**
 * @file MpiConverterBase.cpp
 * @brief Source of the Base class for MPI data converter
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/MpiConverterBase.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   MpiConverterBase::MpiConverterBase()
      : IConverter(), mNeedEmptyComm(false), mPacks(0)
   {
   }

   MpiConverterBase::~MpiConverterBase()
   {
      // Cleanup mpi datatypes
      void cleanupTypes();
   }

   void MpiConverterBase::resetFwdPositions()
   {
      this->mspFBuffers->resetPositions();
   }

   void MpiConverterBase::resetBwdPositions()
   {
      this->mspBBuffers->resetPositions();
   }

   void MpiConverterBase::setBuffers(typename MpiConverterBase::SharedFwdBufferType spFwd, typename MpiConverterBase::SharedBwdBufferType spBwd)
   {
      // Set the forward buffers
      this->mspFBuffers = spFwd;

      // Set the backward buffers
      this->mspBBuffers = spBwd;
   }

   void MpiConverterBase::cleanupTypes()
   {
      #if defined QUICC_MPIPACK_MANUAL
         // No cleanup is required

      #else
         // Cleanup the Fwd Types
         for(auto it = this->mFTypes.begin(); it != this->mFTypes.end(); ++it)
         {
            MPI_Type_free(&(*it));
         }

         // Cleanup the Bwd Types
         for(auto it = this->mBTypes.begin(); it != this->mBTypes.end(); ++it)
         {
            MPI_Type_free(&(*it));
         }
      #endif //defined QUICC_MPIPACK_MANUAL
   }

   void MpiConverterBase::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat memComm = 0.0;

      // General communication storage
      memComm += Debug::MemorySize<int>::BYTES*(this->mFCpuGroup.size());
      memComm += Debug::MemorySize<int>::BYTES*(this->mBCpuGroup.size());
      memComm += Debug::MemorySize<int>::BYTES*(this->mFSizes.size());
      memComm += Debug::MemorySize<int>::BYTES*(this->mBSizes.size());
      memComm += Debug::MemorySize<int>::BYTES*(this->mForwardPacks.size());
      memComm += Debug::MemorySize<int>::BYTES*(this->mBackwardPacks.size());

      StorageProfilerMacro_update(StorageProfilerMacro::MPI, memComm);
      StorageProfilerMacro_update(StorageProfilerMacro::MPICOMM, memComm);

      MHDFloat memTypes = 0.0;

      // MPI datatypes storage
      memTypes += Debug::MemorySize<int>::BYTES*(this->mFTypes.size());
      memTypes += Debug::MemorySize<int>::BYTES*(this->mBTypes.size());

      StorageProfilerMacro_update(StorageProfilerMacro::MPI, 22*memTypes);
      StorageProfilerMacro_update(StorageProfilerMacro::MPITYPES, 22*memTypes);
#endif // QUICC_STORAGEPROFILE
   }

}
}
