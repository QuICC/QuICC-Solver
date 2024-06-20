/**
 * @file MpiConverterSendRecv.cpp
 * @brief Source of the Send/Receive MPI data converter
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/MpiConverterSendRecv.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"

namespace QuICC {

namespace Parallel {

   MpiConverterSendRecv::MpiConverterSendRecv()
      : mIsSending(false), mIsReceiving(false), mActiveSend(0), mActiveReceive(0)
   {
   }

   MpiConverterSendRecv::~MpiConverterSendRecv()
   {
      // Cleanup the requests memory
      this->cleanupRequests();
   }

   void MpiConverterSendRecv::setup()
   {
      // Initialize positions
      this->resetFwdPositions();
      this->resetBwdPositions();

      // Initialise the active packs
      this->mActiveSend = 0;
      this->mActiveReceive = 0;
      this->mActiveDirection = TransformDirection::FORWARD;

      // setup the communication requests
      this->setupRequests();
   }

   void MpiConverterSendRecv::setupCommunication(const int packs, const TransformDirection::Id direction)
   {
      // Store the number of packs in active transfer
      this->mActiveSend = this->mPacks;
      this->mActiveReceive = this->mPacks;
      this->mActiveDirection = this->mDirection;

      // Store the number of packs in the next communication in direction
      this->mPacks = packs;
      this->mDirection = direction;
   }

   void MpiConverterSendRecv::prepareForwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Make sure the buffer is free
         this->syncFwdBuffer();

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nFCpu(), this->pRecvFRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         QuICCEnv().check(ierr, 712);

         // Can only happen if something went really wrong
         if(!flag)
         {
            QuICCEnv().abort("Setting up MPI requests for forward received failed for ID = " + std::to_string(static_cast<int>(this->mTraId)));
         }

         // Prepost the receive calls
         ierr = MPI_Startall(this->nFCpu(), this->pRecvFRequests(this->mPacks));
         QuICCEnv().check(ierr, 722);
         this->resetBwdPositions();
         this->mIsReceiving = true;
      }
   }

   void MpiConverterSendRecv::initiateBackwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
      	// Synchronize
         QuICCEnv().synchronize(this->mTraId);

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nBCpu(), this->pSendBRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         QuICCEnv().check(ierr, 713);

         // Can only happen if something went really wrong
         if(!flag)
         {
            QuICCEnv().abort("Starting MPI communicationr for forward send requests failed for ID = " + std::to_string(static_cast<int>(this->mTraId)));
         }

         // Post non blocking send calls
         ierr = MPI_Startall(this->nBCpu(), this->pSendBRequests(this->mPacks));
         QuICCEnv().check(ierr, 723);
         this->resetFwdPositions();
         this->mIsSending = true;
      }
   }

   void MpiConverterSendRecv::prepareBackwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Make sure the buffer is free
         this->syncBwdBuffer();

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nBCpu(), this->pRecvBRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         QuICCEnv().check(ierr, 714);

         // Can only happen if something went really wrong
         if(!flag)
         {
            QuICCEnv().abort("Setting up MPI requests for backward received failed for ID = " + std::to_string(static_cast<int>(this->mTraId)));
         }

         // Prepost the receive calls
         ierr = MPI_Startall(this->nBCpu(), this->pRecvBRequests(this->mPacks));
         QuICCEnv().check(ierr, 724);
         this->resetFwdPositions();
         this->mIsReceiving = true;
      }
   }

   void MpiConverterSendRecv::initiateForwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
      	// Synchronize
         QuICCEnv().synchronize(mTraId);

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nFCpu(), this->pSendFRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         QuICCEnv().check(ierr, 715);

         // Can only happen if something went really wrong
         if(!flag)
         {
            QuICCEnv().abort("Starting MPI communicationr for forward send requests failed for ID = " + std::to_string(static_cast<int>(this->mTraId)));
         }

         // Post non blocking send calls
         ierr = MPI_Startall(this->nFCpu(), this->pSendFRequests(this->mPacks));
         QuICCEnv().check(ierr, 725);
         this->resetBwdPositions();
         this->mIsSending = true;
      }
   }

   void MpiConverterSendRecv::syncFwdBuffer()
   {
      if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == TransformDirection::BACKWARD)
      {
         this->syncFwdRecvBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == TransformDirection::FORWARD)
      {
         this->syncFwdSendBuffer();

      } else if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncFwdRecvBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncFwdSendBuffer();
      }
   }

   void MpiConverterSendRecv::syncBwdBuffer()
   {
      if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == TransformDirection::BACKWARD)
      {
         this->syncBwdSendBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == TransformDirection::FORWARD)
      {
         this->syncBwdRecvBuffer();

      } else if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncBwdSendBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncBwdRecvBuffer();
      }
   }

   void MpiConverterSendRecv::syncFwdRecvBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveReceive > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nFCpu(), this->pRecvFRequests(this->mActiveReceive), &flag, MPI_STATUSES_IGNORE);
         QuICCEnv().check(ierr, 731);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nFCpu(), this->pRecvFRequests(this->mActiveReceive), MPI_STATUSES_IGNORE);
            QuICCEnv().check(ierr, 741);
         }

         // Clear active packs
         this->mActiveReceive = 0;
      }
   }

   void MpiConverterSendRecv::syncBwdRecvBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveReceive > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nBCpu(), this->pRecvBRequests(this->mActiveReceive), &flag, MPI_STATUSES_IGNORE);
        	QuICCEnv().check(ierr, 732);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nBCpu(), this->pRecvBRequests(this->mActiveReceive), MPI_STATUSES_IGNORE);
        	   QuICCEnv().check(ierr, 742);
         }

         // Clear active packs
         this->mActiveReceive = 0;
      }
   }

   void MpiConverterSendRecv::syncFwdSendBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveSend > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nBCpu(), this->pSendFRequests(this->mActiveSend), &flag, MPI_STATUSES_IGNORE);
        	QuICCEnv().check(ierr, 733);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nBCpu(), this->pSendFRequests(this->mActiveSend), MPI_STATUSES_IGNORE);
        	   QuICCEnv().check(ierr, 743);
         }

         // Clear active packs
         this->mActiveSend = 0;
      }
   }

   void MpiConverterSendRecv::syncBwdSendBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveSend > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nFCpu(), this->pSendBRequests(this->mActiveSend), &flag, MPI_STATUSES_IGNORE);
        	QuICCEnv().check(ierr, 734);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nFCpu(), this->pSendBRequests(this->mActiveSend), MPI_STATUSES_IGNORE);
        	   QuICCEnv().check(ierr, 744);
         }

         // Clear active packs
         this->mActiveSend = 0;
      }
   }

   int MpiConverterSendRecv::sendDest(const int id, const int ref, const int size) const
   {
      // Create send ring
      return ((id + 1 + ref) % size);
   }

   int MpiConverterSendRecv::recvSrc(const int id, const int ref, const int size) const
   {
      // Create recv ring
      return ((size - 1 - id + ref) % size);
   }

   void MpiConverterSendRecv::setupRequests()
   {
      // Storage for global location flags
      int dest;
      int src;
      int tag;
      // Storage for CPU group location flags
      int grpMe;
      int grpDest;
      int grpSrc;

      // Shift tag to produce unique tag in 2D distribution
      // (Probably not required anymore but doesn't harm)
      int tagShift = 0;
      if(this->mTraId == Dimensions::Transform::TRA1D)
      {
         tagShift = 0;
      } else if(this->mTraId == Dimensions::Transform::TRA2D)
      {
         tagShift = QuICCEnv().size();
      } else
      {
         QuICCEnv().abort("Settig up MPI requests failed");
      }

      // MPI error code
      int ierr;

      // Storage for the number of packs
      int packs;

      // Initialise forward transform requests
      for(int k = 0; k < this->mForwardPacks.size(); ++k)
      {
         // Get the pack size
         packs = this->mForwardPacks(k);

         // Initialise receive forward with empty requests
         this->mRecvFRequests.insert(std::make_pair(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            this->mRecvFRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Initialise send backward with empty requests
         this->mSendBRequests.insert(std::make_pair(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            this->mSendBRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Create receive forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), QuICCEnv().id(this->mTraId)));

            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nFCpu());

            // Get source MPI rank in group
            src = this->fCpu(grpSrc);

            // Set shifted MPI tag to make it unique
            tag = src + tagShift;

            //Safety asserts
            assert(static_cast<size_t>(grpSrc) < this->mFSizes.size());
            assert(static_cast<size_t>(grpSrc) < this->mRecvFRequests.at(packs).size());

            // initialise the Recv request
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Recv_init(this->mspFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MpiTypes::type<FwdBufferType::DataType>(), src, tag, QuICCEnv().comm(this->mTraId), &(this->mRecvFRequests.at(packs).at(grpSrc)));
               QuICCEnv().check(ierr, 981);
            #else
               ierr = MPI_Recv_init(this->mspFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MPI_PACKED, src, tag, QuICCEnv().comm(this->mTraId), &(this->mRecvFRequests.at(packs).at(grpSrc)));
               QuICCEnv().check(ierr, 981);
            #endif //defined QUICC_MPIPACK_MANUAL
         }

         // Create send backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), QuICCEnv().id(this->mTraId)));

            // Set shifted MPI tag to make it unique
            tag = QuICCEnv().id(this->mTraId) + tagShift;

            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nBCpu());

            // Get destination MPI rank in group
            dest = this->bCpu(grpDest);

            //Safety asserts
            assert(static_cast<size_t>(grpDest) < this->mBSizes.size());
            assert(static_cast<size_t>(grpDest) < this->mSendBRequests.at(packs).size());

            // initialise the Send request
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Send_init(this->mspBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MpiTypes::type<BwdBufferType::DataType>(), dest, tag, QuICCEnv().comm(this->mTraId), &(this->mSendBRequests.at(packs).at(grpDest)));
               QuICCEnv().check(ierr, 982);
            #else
               ierr = MPI_Send_init(this->mspBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MPI_PACKED, dest, tag, QuICCEnv().comm(this->mTraId), &(this->mSendBRequests.at(packs).at(grpDest)));
               QuICCEnv().check(ierr, 982);
            #endif //defined QUICC_MPIPACK_MANUAL
         }
      }

      // Initialise backward transform requests
      for(int k = 0; k < this->mBackwardPacks.size(); ++k)
      {
         // Get the pack size
         packs = this->mBackwardPacks(k);

         // Initialise receive backward with empty requests
         this->mRecvBRequests.insert(std::make_pair(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            this->mRecvBRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Initialise send forward with empty requests
         this->mSendFRequests.insert(std::make_pair(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            this->mSendFRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Create receive backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), QuICCEnv().id(this->mTraId)));

            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nBCpu());

            // Get source MPI rank in group
            src = this->bCpu(grpSrc);

            // Set shifted MPI tag to make it unique
            tag = src + tagShift;

            // initialise the Recv request
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Recv_init(this->mspBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MpiTypes::BwdBufferType::DataType>(), src, tag, QuICCEnv().comm(this->mTraId), &(this->mRecvBRequests.at(packs).at(grpSrc)));
            #else
               ierr = MPI_Recv_init(this->mspBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MPI_PACKED, src, tag, QuICCEnv().comm(this->mTraId), &(this->mRecvBRequests.at(packs).at(grpSrc)));
            #endif //defined QUICC_MPIPACK_MANUAL
            QuICCEnv().check(ierr, 983);
         }

         // Create send forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), QuICCEnv().id(this->mTraId)));

            // Set shifted MPI tag to make it unique
            tag = QuICCEnv().id(this->mTraId) + tagShift;

            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nFCpu());

            // Get destination MPI rank in group
            dest = this->fCpu(grpDest);

            // initialise the Send request
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Send_init(this->mspFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MpiTypes::type<FwdBufferType::DataType>(), dest, tag, QuICCEnv().comm(this->mTraId), &(this->mSendFRequests.at(packs).at(grpDest)));
            #else
               ierr = MPI_Send_init(this->mspFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MPI_PACKED, dest, tag, QuICCEnv().comm(this->mTraId), &(this->mSendFRequests.at(packs).at(grpDest)));
               QuICCEnv().check(ierr, 984);
            #endif //defined QUICC_MPIPACK_MANUAL
         }
      }
   }

   void MpiConverterSendRecv::cleanupRequests()
   {
      // Free requests from Recv B
      for(auto it = this->mRecvBRequests.begin(); it != this->mRecvBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Recv F
      for(auto it = this->mRecvFRequests.begin(); it != this->mRecvFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send B
      for(auto it = this->mSendBRequests.begin(); it != this->mSendBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send F
      for(auto it = this->mSendFRequests.begin(); it != this->mSendFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }
   }

   void MpiConverterSendRecv::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MpiConverterBase::profileStorage();

      MHDFloat memComm = 0.0;

      // Requests communication storage
      memComm += Debug::MemorySize<int>::BYTES*5.0;
      MHDFloat memReq = (Debug::MemorySize<int>::BYTES + Debug::MemorySize<int>::BYTES);
      memComm += memReq*(this->mRecvBRequests.size());
      memComm += memReq*(this->mRecvFRequests.size());
      memComm += memReq*(this->mSendBRequests.size());
      memComm += memReq*(this->mSendFRequests.size());

      StorageProfilerMacro_update(StorageProfilerMacro::MPICOMM, memComm);
      StorageProfilerMacro_update(StorageProfilerMacro::MPI, memComm);
#endif // QUICC_STORAGEPROFILE
   }

}
}
