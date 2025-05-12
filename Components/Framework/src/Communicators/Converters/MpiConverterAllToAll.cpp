/**
 * @file MpiConverter.cpp
 * @brief Source of the AlltoAll MPI converter
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/MpiConverterAllToAll.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   MpiConverterAllToAll::MpiConverterAllToAll()
   {
      // Keep empty communications
      this->mNeedEmptyComm = true;
   }

   MpiConverterAllToAll::~MpiConverterAllToAll()
   {
   }

   void MpiConverterAllToAll::setup()
   {
      // Initialize positions
      this->resetFwdPositions();
      this->resetBwdPositions();

      // setup the communication data
      this->setupCommData();
   }

   void MpiConverterAllToAll::setupCommunication(const int packs, const TransformDirection::Id direction)
   {
      // Store the number of packs in the next communication in direction
      this->mPacks = packs;
      this->mDirection = direction;
   }

   void MpiConverterAllToAll::prepareForwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         this->resetBwdPositions();
      }
   }

   void MpiConverterAllToAll::initiateBackwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         this->resetFwdPositions();
         this->setFwdCommSizes();
         this->setBwdCommSizes();

         // All-to-all data exchange
         #if defined QUICC_MPIPACK_MANUAL
            MPI_Alltoallv(this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), this->mBCommTypes[0], this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), this->mFCommTypes[0], QuICCEnv().comm(this->mTraId));
         #else
            MPI_Alltoallw(this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), &this->mBCommTypes[0], this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), &this->mFCommTypes[0], QuICCEnv().comm(this->mTraId));
         #endif //defined QUICC_MPIPACK_MANUAL
      }
   }

   void MpiConverterAllToAll::prepareBackwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         this->resetFwdPositions();
      }
   }

   void MpiConverterAllToAll::initiateForwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         this->resetBwdPositions();
         this->setFwdCommSizes();
         this->setBwdCommSizes();

         // All-to-all data exchange
         #if defined QUICC_MPIPACK_MANUAL
            MPI_Alltoallv(this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), this->mFCommTypes[0], this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), this->mBCommTypes[0], QuICCEnv().comm(this->mTraId));
         #else
            MPI_Alltoallw(this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), &this->mFCommTypes[0], this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), &this->mBCommTypes[0], QuICCEnv().comm(this->mTraId));
         #endif //defined QUICC_MPIPACK_MANUAL
      }
   }

   void MpiConverterAllToAll::setFwdCommSizes()
   {
      // All-to-all data FWD sizes
      for(size_t i = 0; i < this->mFSizes.size(); i++)
      {
         this->mFCommSizes(i) = this->mPacks*this->mFSizes.at(i);
      }
   }

   void MpiConverterAllToAll::setBwdCommSizes()
   {
      // All-to-all data FWD sizes
      for(size_t i = 0; i < this->mBSizes.size(); i++)
      {
         this->mBCommSizes(i) = this->mPacks*this->mBSizes.at(i);
      }
   }

   void MpiConverterAllToAll::setupCommData()
   {
      #if defined QUICC_MPIPACK_MANUAL
         this->mFCommTypes.push_back(MpiTypes::type<FwdBufferType::DataType>());

         this->mBCommTypes.push_back(MpiTypes::type<BwdBufferType::DataType>());
      #else
         for(size_t i = 0; i < this->mFSizes.size(); i++)
         {
            this->mFCommTypes.push_back(MPI_PACKED);
         }

         for(size_t i = 0; i < this->mBSizes.size(); i++)
         {
            this->mBCommTypes.push_back(MPI_PACKED);
         }
      #endif //defined QUICC_MPIPACK_MANUAL

      this->mFCommSizes.resize(this->mFSizes.size());
      this->mBCommSizes.resize(this->mBSizes.size());
   }

   void MpiConverterAllToAll::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MpiConverterBase::profileStorage();

      MHDFloat memComm = 0.0;

      // Requests communication storage
      memComm += Debug::MemorySize<int>::BYTES*(this->mFCommSizes.size());
      memComm += Debug::MemorySize<int>::BYTES*(this->mBCommSizes.size());

      StorageProfilerMacro_update(StorageProfilerMacro::MPI, memComm);
      StorageProfilerMacro_update(StorageProfilerMacro::MPICOMM, memComm);

      MHDFloat memTypes = 0.0;

      // MPI datatypes storage
      memTypes += Debug::MemorySize<int>::BYTES*(this->mFCommTypes.size());
      memTypes += Debug::MemorySize<int>::BYTES*(this->mBCommTypes.size());

      StorageProfilerMacro_update(StorageProfilerMacro::MPI, 22*memTypes);
      StorageProfilerMacro_update(StorageProfilerMacro::MPITYPES, 22*memTypes);
#endif // QUICC_STORAGEPROFILE
   }

}
}
