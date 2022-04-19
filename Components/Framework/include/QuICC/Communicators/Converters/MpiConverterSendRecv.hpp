/**
 * @file MpiConverterSendRecv.hpp
 * @brief Implementation of the MPI data converter with send/recv communication
 */

#ifndef QUICC_PARALLEL_MPICONVERTERSENDRECV_HPP
#define QUICC_PARALLEL_MPICONVERTERSENDRECV_HPP

// Configuration includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// System includes
//
#include <type_traits>

// External includes
//

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/Communicators/Converters/MpiConverterBase.hpp"
#include "QuICC/Communicators/Converters/MpiConverterTools.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the MPI data converter with send/recv communication.
    */
   class MpiConverterSendRecv: public MpiConverterBase
   {
      public:
         /// Typedef for forward datatype
         typedef MpiConverterBase::RealFwdData RealFwdData;

         /// Typedef for forward datatype
         typedef MpiConverterBase::ComplexFwdData ComplexFwdData;

         /// Typedef for backward datatype
         typedef MpiConverterBase::RealBwdData RealBwdData;

         /// Typedef for backward datatype
         typedef MpiConverterBase::ComplexBwdData ComplexBwdData;

         /**
          * @brief Constructor
          */
         MpiConverterSendRecv();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterSendRecv();

         /**
          * @brief Finish the setup of the converter
          */
         virtual void setup() override;

         /**
          * @brief Setup upcoming communication
          *
          * @param packs Number of packets in communication packing
          */
         virtual void setupCommunication(const int packs, const TransformDirection::Id direction) override;

         /**
          * @brief Start persistent send for forward transform
          */
         virtual void initiateForwardSend() override;

         /**
          * @brief Post persistent receive for forward transform
          */
         virtual void prepareForwardReceive() override;

         /**
          * @brief Start persistent send for backward transform
          */
         virtual void initiateBackwardSend() override;

         /**
          * @brief Post persistent receive for backward transform
          */
         virtual void prepareBackwardReceive() override;

         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const override;

      protected:
         /**
          * @brief Send forward data
          *
          * @param data Data to send
          */
         template <typename TFwd> void sendFwd(const TFwd &data);

         /**
          * @brief Send backward data
          *
          * @param data Data to send
          */
         template <typename TBwd> void sendBwd(const TBwd &data);

         /**
          * @brief Receive forward data
          *
          * @param rData Storage for received data
          */
         template <typename TFwd> void receiveFwd(TFwd &rData);

         /**
          * @brief Receive backward data
          *
          * @param rData Storage for received data
          */
         template <typename TBwd> void receiveBwd(TBwd &rData);

      private:
         /**
          * @brief Get a pointer to the receive backward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pRecvBRequests(const int size);

         /**
          * @brief Get a pointer to the receive forward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pRecvFRequests(const int size);

         /**
          * @brief Get a pointer to the send backward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pSendBRequests(const int size);

         /**
          * @brief Get a pointer to the send forward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pSendFRequests(const int size);

         /**
          * @brief Make sure forward buffer is available
          */
         void syncFwdBuffer();

         /**
          * @brief Make sure backward buffer is available
          */
         void syncBwdBuffer();

         /**
          * @brief Make sure forward buffer used by send is available
          */
         void syncFwdSendBuffer();

         /**
          * @brief Make sure backward buffer used by send is available
          */
         void syncBwdSendBuffer();

         /**
          * @brief Make sure forward buffer used by receive is available
          */
         void syncFwdRecvBuffer();

         /**
          * @brief Make sure backward buffer used by receive is available
          */
         void syncBwdRecvBuffer();

         /**
          * @brief Get ring recv source for id
          *
          * @param id   ID of the CPU
          * @param ref  ID of the reference CPU
          * @param size Size of the CPU group
          */
         int  recvSrc(const int id, const int ref, const int size) const;

         /**
          * @brief Get ring send destination for id
          *
          * @param id   ID of the CPU
          * @param ref  ID of the reference CPU
          * @param size Size of the CPU group
          */
         int  sendDest(const int id, const int ref, const int size) const;

         /**
          * @brief Setup the MPI communication requests requests
          */
         void setupRequests();

         /**
          * @brief Cleanup the MPI communication requests
          */
         void cleanupRequests();

         /**
          * @brief Sending communication status
          */
         bool  mIsSending;

         /**
          * @brief Receiving communication status
          */
         bool  mIsReceiving;

         /**
          * @brief Direction of active operation
          */
         TransformDirection::Id   mActiveDirection;

         /**
          * @brief The number of packs in the "previous/active" send operations
          */
         int   mActiveSend;

         /**
          * @brief The number of packs in the "previous/active" receive operations
          */
         int   mActiveReceive;

         /**
          * @brief Storage for the non blocking communication requests: Recv F
          */
         std::map<int, std::vector<MPI_Request> >  mRecvFRequests;

         /**
          * @brief Storage for the non blocking communication requests: Recv B
          */
         std::map<int, std::vector<MPI_Request> >  mRecvBRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send F
          */
         std::map<int, std::vector<MPI_Request> >  mSendFRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send B
          */
         std::map<int, std::vector<MPI_Request> >  mSendBRequests;
   };

   inline MPI_Request * MpiConverterSendRecv::pRecvBRequests(const int size)
   {
      return &(this->mRecvBRequests.at(size).front());
   }

   inline MPI_Request * MpiConverterSendRecv::pRecvFRequests(const int size)
   {
      return &(this->mRecvFRequests.at(size).front());
   }

   inline MPI_Request * MpiConverterSendRecv::pSendBRequests(const int size)
   {
      return &(this->mSendBRequests.at(size).front());
   }

   inline MPI_Request * MpiConverterSendRecv::pSendFRequests(const int size)
   {
      return &(this->mSendFRequests.at(size).front());
   }

   template <typename TFwd> void MpiConverterSendRecv::sendFwd(const TFwd& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

      // Make sure the buffer is free
      this->syncFwdBuffer();

      // Pack data into send buffer
      for(int id = 0; id < this->nFCpu(); ++id)
      {
         #if defined QUICC_MPIPACK_MANUAL
            MpiConverterTools::pack(this->mspFBuffers->at(id), this->mspFBuffers->pos(id), data, this->mFTypes.at(id));
         #else
            int ierr = MPI_Pack(const_cast<typename TFwd::PointType *>(data.data().data()), 1, this->mFTypes.at(id), this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mspFBuffers->pos(id)), MpiFramework::transformComm(this->mTraId));
            QuICCEnv().check(ierr, 761);
         #endif //defined QUICC_MPIPACK_MANUAL
      }
   }

   template <typename TBwd> void MpiConverterSendRecv::sendBwd(const TBwd& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

      // Make sure the buffer is free
      this->syncBwdBuffer();

      // Pack data into send buffer
      for(int id = 0; id < this->nBCpu(); ++id)
      {
         #if defined QUICC_MPIPACK_MANUAL
            MpiConverterTools::pack(this->mspBBuffers->at(id), this->mspBBuffers->pos(id), data, this->mBTypes.at(id));
         #else
            int ierr = MPI_Pack(const_cast<typename TBwd::PointType *>(data.data().data()), 1, this->mBTypes.at(id), this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mspBBuffers->pos(id)), MpiFramework::transformComm(this->mTraId));
            QuICCEnv().check(ierr, 762);
         #endif //defined QUICC_MPIPACK_MANUAL
      }
   }

   template <typename TFwd> void MpiConverterSendRecv::receiveFwd(TFwd &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

      // Communication interface is receiving the data
      if(this->mIsReceiving)
      {
         // Number of receive calls in total
         int keepWaiting = this->nFCpu();
         int count = 0;
         ArrayI   idx = ArrayI::Zero(this->nFCpu());

         // Wait until everything has been received
         while(keepWaiting != 0)
         {
            // Wait for some of the requests to finish
            #ifdef QUICC_DEBUG
               MPI_Status stats[this->nFCpu()];
               int ierr = MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), stats);
               QuICCEnv().check(ierr, 771);
               DebuggerMacro_msg("Received FWD packs", 5);
            #else
               int ierr = MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);
               QuICCEnv().check(ierr, 771);
            #endif //QUICC_DEBUG

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               DebuggerMacro_showValue("Tag: ", 6, stats[id].MPI_TAG);
               DebuggerMacro_showValue("-> From: ", 6, stats[id].MPI_SOURCE);
               int pos = idx(id);

               #if defined QUICC_MPIPACK_MANUAL
                  MpiConverterTools::unpack(rData, this->mFTypes.at(pos), this->mspFBuffers->at(pos), this->mspFBuffers->pos(pos));
               #else
                  int ierr = MPI_Unpack(this->mspFBuffers->at(pos), this->sizeFPacket(pos), &(this->mspFBuffers->pos(pos)), rData.rData().data(), 1, this->mFTypes.at(pos), MpiFramework::transformComm(this->mTraId));
                  QuICCEnv().check(ierr, 763);
               #endif //defined QUICC_MPIPACK_MANUAL
            }

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;

      // Data is here and just need to be unpacked
      } else
      {
         // Unpack data from receive buffer
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            DebuggerMacro_msg("Unpacking FWD packs", 5);

            #if defined QUICC_MPIPACK_MANUAL
               MpiConverterTools::unpack(rData, this->mFTypes.at(id), this->mspFBuffers->at(id), this->mspFBuffers->pos(id));
            #else
               int ierr = MPI_Unpack(this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mspFBuffers->pos(id)), rData.rData().data(), 1, this->mFTypes.at(id), MpiFramework::transformComm(this->mTraId));
               QuICCEnv().check(ierr, 764);
            #endif //defined QUICC_MPIPACK_MANUAL

         }
      }
   }

   template <typename TBwd> void MpiConverterSendRecv::receiveBwd(TBwd &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

      // Communication interface is receiving the data
      if(this->mIsReceiving)
      {
         // Number of receive calls in total
         int keepWaiting = this->nBCpu();
         int count = 0;
         ArrayI   idx = ArrayI::Zero(this->nBCpu());

         // Wait until everything has been received
         while(keepWaiting != 0)
         {
            // Wait for some of the requests to finish
            #ifdef QUICC_DEBUG
               MPI_Status stats[this->nBCpu()];
               int ierr = MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), stats);
               QuICCEnv().check(ierr, 772);
               DebuggerMacro_msg("Received BWD packs", 5);
            #else
               int ierr = MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);
               QuICCEnv().check(ierr, 772);
            #endif //QUICC_DEBUG

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               DebuggerMacro_showValue("Tag: ", 6, stats[id].MPI_TAG);
               DebuggerMacro_showValue("-> From: ", 6, stats[id].MPI_SOURCE);
               int pos = idx(id);

               #if defined QUICC_MPIPACK_MANUAL
                  MpiConverterTools::unpack(rData, this->mBTypes.at(pos), this->mspBBuffers->at(pos), this->mspBBuffers->pos(pos));
               #else
                  int ierr = MPI_Unpack(this->mspBBuffers->at(pos), this->sizeBPacket(pos), &(this->mspBBuffers->pos(pos)), rData.rData().data(), 1, this->mBTypes.at(pos), MpiFramework::transformComm(this->mTraId));
                  QuICCEnv().check(ierr, 765);
               #endif //defined QUICC_MPIPACK_MANUAL
            }

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;

      // Data is here and just need to be unpacked
      } else
      {
         // Unpack data from receive buffer
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            DebuggerMacro_msg("Unpacking BWD packs", 5);

            #if defined QUICC_MPIPACK_MANUAL
               MpiConverterTools::unpack(rData, this->mBTypes.at(id), this->mspBBuffers->at(id), this->mspBBuffers->pos(id));
            #else
               int ierr = MPI_Unpack(this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mspBBuffers->pos(id)), rData.rData().data(), 1, this->mBTypes.at(id), MpiFramework::transformComm(this->mTraId));
               QuICCEnv().check(ierr, 766);
            #endif //defined QUICC_MPIPACK_MANUAL
         }
      }
   }

}
}

#endif // QUICC_PARALLEL_MPICONVERTERSENDRECV_HPP
