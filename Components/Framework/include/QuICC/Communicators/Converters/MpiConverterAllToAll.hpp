/**
 * @file MpiConverterAllToAll.hpp
 * @brief Implementation of the MPI data converter with AllToAll data exchange
 */

#ifndef QUICC_PARALLEL_MPICONVERTERALLTOALL_HPP
#define QUICC_PARALLEL_MPICONVERTERALLTOALL_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/Communicators/Converters/MpiConverterBase.hpp"
#include "QuICC/Communicators/Converters/MpiConverterTools.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the MPI data converter with AllToAll data exchange.
    */
   class MpiConverterAllToAll: public MpiConverterBase
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
         MpiConverterAllToAll();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterAllToAll();

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
          * @brief Set forward communication sizes
          */
         void setFwdCommSizes();

         /**
          * @brief Set backward communication sizes
          */
         void setBwdCommSizes();

         /**
          * @brief Setup the communication data
          */
         void setupCommData();

         /**
          * @brief Forward communication sizes
          */
         ArrayI mFCommSizes;

         /**
          * @brief backward communication sizes
          */
         ArrayI mBCommSizes;

         /**
          * @brief MPI datatypes for forward communication
          */
         std::vector<MPI_Datatype>  mFCommTypes;

         /**
          * @brief MPI datatypes for backward communication
          */
         std::vector<MPI_Datatype>  mBCommTypes;
   };

   template <typename TFwd> void MpiConverterAllToAll::sendFwd(const TFwd& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

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

   template <typename TBwd> void MpiConverterAllToAll::sendBwd(const TBwd& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

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

   template <typename TFwd> void MpiConverterAllToAll::receiveFwd(TFwd &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

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

   template <typename TBwd> void MpiConverterAllToAll::receiveBwd(TBwd &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

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

#endif // QUICC_PARALLEL_MPICONVERTERALLTOALL_HPP
