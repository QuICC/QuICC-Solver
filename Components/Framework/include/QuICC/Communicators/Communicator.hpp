/**
 * @file Communicator.hpp
 * @brief Implementation of a 1D communicator
 */

#ifndef QUICC_PARALLEL_COMMUNICATOR_HPP
#define QUICC_PARALLEL_COMMUNICATOR_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Arithmetics/Set.hpp"
#include "QuICC/Arithmetics/SetNeg.hpp"
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Communicators/CommunicationBuffer.hpp"
#include "QuICC/Communicators/CommunicatorStorage.hpp"
#include "QuICC/Communicators/Converters/SerialConverter.hpp"
#ifdef QUICC_MPI
#include "QuICC/Communicators/Converters/MpiConverter.hpp"
#endif // QUICC_MPI

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of 1D communicator
    */
   class Communicator: public CommunicatorStorage
   {
      public:
         /**
         * @brief Constructor
         */
         Communicator();

         /**
         * @brief Destructor
         */
         ~Communicator();

         /**
          * @brief Initialise the communicator for given dimension
          *
          * @param spSetupFwd Setup object for the forward 1D type
          * @param spSetupBwd Setup object for the backward 1D type
          */
         void init(const Dimensions::Transform::Id id, SharedFwdSetupType spSetupFwd, SharedBwdSetupType spSetupBwd);

         /**
          * @brief Initialise the 2D/3D converter
          *
          * @param spRes      Shared resolution information
          * @param packs1DFwd Packs information for first forward exchange
          * @param packs1DBwd Packs information for first backward exchange
          * @param packs2DFwd Packs information for second forward exchange
          * @param packs2DBwd Packs information for second backward exchange
          * @param split      Location where the MPI splitting takes place
          */
         void initConverter(SharedResolution spRes, const std::vector<ArrayI>& packs, Splitting::Locations::Id split);

         /**
          * @brief Transfer forward data to next step
          */
         template <typename T> void transferForward(Dimensions::Transform::Id tId, T& rData, const bool freeData = true);

         /**
          * @brief Transfer backward data to next step
          */
         template <typename T> void transferBackward(Dimensions::Transform::Id tId, T& rData);

         /**
          * @brief Receive forward the data
          */
         template <typename T> void receiveForward(Dimensions::Transform::Id tId, T& pData);

         /**
          * @brief Receive backward the data
          */
         template <typename T> void receiveBackward(Dimensions::Transform::Id tId, T& pData);

         /**
          * @brief Provide physical storage
          */
         template <typename T> void providePhysical(T& pData);

         /**
          * @brief Hold final physical storage
          */
         template <typename T> void holdPhysical(T& rData);

      protected:

      private:
         /**
          * @brief Create serial converter
          *
          * @param spRes      Shared resolution information
          */
         template <Dimensions::Transform::Id TId> void createSerialConverter(SharedResolution spRes);

#ifdef QUICC_MPI
         /**
          * @brief Create MPI converter
          *
          * @param spRes      Shared resolution information
          */
         template <Dimensions::Transform::Id TId> void createMpiConverter(SharedResolution spRes, const ArrayI& packsFwd, const ArrayI& packsBwd);
#endif // QUICC_MPI

   };

   template <typename T>
      void Communicator::receiveForward(Dimensions::Transform::Id tId, T& pData)
   {
      // Debugger message
      DebuggerMacro_msg("receiveForward ID = " + std::to_string(static_cast<int>(tId)), 5);

      if(this->hasConverter(Dimensions::jump(tId,1)))
      {
         DebuggerMacro_msg("has converter", 6);
         this->converter(Dimensions::jump(tId,1)).getFwd(pData, this->storage(tId));
      } else
      {
         DebuggerMacro_msg("without converter", 6);
         this->storage(tId).recoverFwd(pData);
      }
   }

   template <typename T>
      void Communicator::receiveBackward(Dimensions::Transform::Id tId, T& pData)
   {
      // Debugger message
      DebuggerMacro_msg("receiveBackward ID = " + std::to_string(static_cast<int>(tId)), 5);

      if(this->hasConverter(tId))
      {
         DebuggerMacro_msg("has converter", 6);
         this->converter(tId).getBwd(pData, this->storage(tId));
      }
      else
      {
         DebuggerMacro_msg("without converter", 6);
         this->storage(tId).recoverBwd(pData);
      }
   }

   template <typename T>
      void Communicator::transferForward(Dimensions::Transform::Id tId, T& rData, const bool freeInput)
   {
      // Debugger message
      DebuggerMacro_msg("transferForward ID = " + std::to_string(static_cast<int>(tId)), 5);

      auto traId = Dimensions::jump(tId,1);
      if(this->hasConverter(traId))
      {
         DebuggerMacro_msg("has converter", 6);
         // Convert data
         this->converter(traId).convertFwd(rData, this->storage(traId));

         // Free input data
         if(freeInput)
         {
            this->storage(tId).freeFwd(rData);
         }
      }
      else
      {
         DebuggerMacro_msg("NEEDS FIX: transferForward is not consistent here", 5);
         DebuggerMacro_msg("transferForward ID = " + std::to_string(static_cast<int>(tId)), 5);
      }
   }

   template <typename T>
      void Communicator::transferBackward(Dimensions::Transform::Id tId, T& rData)
   {
      // Debugger message
      DebuggerMacro_msg("transferBackward ID = " + std::to_string(static_cast<int>(tId)), 5);

      if(this->hasConverter(tId))
      {
         DebuggerMacro_msg("has converter", 6);
         // Convert data
         this->converter(tId).convertBwd(rData, this->storage(Dimensions::jump(tId,-1)));

         // Free input data
         this->storage(tId).freeBwd(rData);
      }
      else
      {
         throw std::logic_error("Something went wrong with the transfer setup");
      }
   }

   template <typename T>
      void Communicator::providePhysical(T& pData)
   {
      this->lastStorage().provideFwd(pData);
   }

   template <typename T>
      void Communicator::holdPhysical(T& rData)
   {
      this->lastStorage().holdFwd(rData);
   }

   template <Dimensions::Transform::Id TId>
      void Communicator::createSerialConverter(SharedResolution spRes)
   {
      // Create shared serial converter
      auto spConv = std::make_shared<SerialConverter>();
      auto spIdxConv = spRes->sim().ss().createIndexConv(TId);

      const Dimensions::Transform::Id jid = Dimensions::jump(TId,-1);
      spIdxConv->init(*spRes, jid);

      // Initialise shared converter
      spConv->init(spRes, jid, spIdxConv);

      // Set 1D/2D converter
      this->addConverter(TId, spConv);
   }

#ifdef QUICC_MPI
      template <Dimensions::Transform::Id TId>
         void Communicator::createMpiConverter(SharedResolution spRes, const ArrayI& packsFwd, const ArrayI& packsBwd)
   {
      // Create shared MPI converter
      auto spConv = std::make_shared<MpiConverter>();
      auto spIdxConv = spRes->sim().ss().createIndexConv(TId);

      const Dimensions::Transform::Id jid = Dimensions::jump(TId,-1);
      spIdxConv->init(*spRes, jid);

      // Initialise the MPI converter
      auto pFTmp = spRes->sim().ss().fwdPtr(jid);
      this-> storage(jid).provideFwd(pFTmp);
      auto pBTmp = spRes->sim().ss().bwdPtr(TId);
      this->storage(TId).provideBwd(pBTmp);
      std::visit(
            [&](auto && pF, auto && pB)
            {
               spConv->init(spRes, jid, *pF, *pB, packsFwd, packsBwd, spIdxConv);
            }
            , pFTmp, pBTmp);
      this->storage(jid).freeFwd(pFTmp);
      this->storage(TId).freeBwd(pBTmp);

      QuICCEnv().synchronize();

      // Create the communication buffers
      auto spBufferOne = std::make_shared<CommunicationBuffer<char> >();
      auto spBufferTwo = std::make_shared<CommunicationBuffer<char> >();

      // Get maximum number of packs
      int pFwd;
      if(packsFwd.size() > 0)
      {
         pFwd = packsFwd(packsFwd.size()-1);
      } else
      {
         pFwd = 0;
      }
      int pBwd;
      if(packsBwd.size() > 0)
      {
         pBwd = packsBwd(packsBwd.size()-1);
      } else
      {
         pBwd = 0;
      }
      int max = std::max(pFwd, pBwd);

      // Allocate first buffers
      spBufferOne->allocate(spConv->fwdSizes(), max);

      // Allocate second buffers
      spBufferTwo->allocate(spConv->bwdSizes(), max);

      // Set communication buffers
      spConv->setBuffers(spBufferOne, spBufferTwo);

      // Set the converter
      this->addConverter(TId, spConv);
   }
#endif // QUICC_MPI

}
}

#endif // QUICC_PARALLEL_COMMUNICATOR_HPP
