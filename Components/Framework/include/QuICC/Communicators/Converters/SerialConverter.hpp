/**
 * @file SerialConverter.hpp
 * @brief Implementation of the serial data converter
 */

#ifndef QUICC_PARALLEL_SERIALCONVERTER_HPP
#define QUICC_PARALLEL_SERIALCONVERTER_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <cassert>
#include <memory>
#include <type_traits>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Enums/TransformDirection.hpp"
#include "QuICC/StorageProviders/DynamicPairProvider.hpp"
#include "QuICC/Communicators/Converters/SerialConverterBase.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the serial data converter.
    */
   class SerialConverter : public SerialConverterBase
   {
      public:
         /// Typedef for forward datatype
         typedef SerialConverterBase::RealFwdData RealFwdData;

         /// Typedef for forward datatype
         typedef SerialConverterBase::ComplexFwdData ComplexFwdData;

         /// Typedef for backward datatype
         typedef SerialConverterBase::RealBwdData RealBwdData;

         /// Typedef for backward datatype
         typedef SerialConverterBase::ComplexBwdData ComplexBwdData;

         /**
          * @brief Constructor
          */
         SerialConverter();

         /**
          * @brief Destructor
          */
         virtual ~SerialConverter();

         /**
          * @brief Initialise the converter
          *
          * @param spRes   Shared Resolution
          * @param id      Dimension index for forward transform
          */
         void init(SharedResolution spRes, const Dimensions::Transform::Id id, std::shared_ptr<IIndexConv> spIdxConv);

         /**
          * @brief Setup the converter
          */
         virtual void setup() override;

         /**
          * @brief Get the converted real data from Bwd to Fwd
          */
         virtual void getFwd(RealFwdData *& pOut, DynamicPairProvider &storage) override;

         /**
          * @brief Get the converted complex data from Bwd to Fwd
          */
         virtual void getFwd(ComplexFwdData *& pOut, DynamicPairProvider &storage) override;

         /**
          * @brief Get the converted real data from Fwd to Bwd
          */
         virtual void getBwd(RealBwdData *& pOut, DynamicPairProvider &storage) override;

         /**
          * @brief Get the converted complex data from Fwd to Bwd
          */
         virtual void getBwd(ComplexBwdData *& pOut, DynamicPairProvider &storage) override;

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

      private:
         /**
          * @brief Get the converted complex data from Bwd to Fwd
          */
         template <typename T> void getFwdImpl(T *& pOut, DynamicPairProvider &storage);

         /**
          * @brief Get the converted real data from Fwd to Bwd
          */
         template <typename T> void getBwdImpl(T *& pOut, DynamicPairProvider &storage);
   };

   template <typename T>
      void SerialConverter::getFwdImpl(T *& pData, DynamicPairProvider  &storage)
   {
      if(pData != nullptr)
      {
         T *pTmp;
         // Recover storage from provider
         storage.recoverFwd(pTmp);

         pData->rData() = pTmp->data();

         storage.freeFwd(*pTmp);
      }
      else
      {
         // Recover storage from provider
         storage.recoverFwd(pData);
      }
   }

   template <typename T>
      void SerialConverter::getBwdImpl(T *& pData, DynamicPairProvider  &storage)
   {
      if(pData != nullptr)
      {
         T *pTmp;

         // Recover storage from provider
         storage.recoverBwd(pTmp);

         pData->rData() = pTmp->data();

         storage.freeBwd(*pTmp);

      } else
      {
         // Recover storage from provider
         storage.recoverBwd(pData);
      }
   }

}
}

#endif // QUICC_PARALLEL_SERIALCONVERTER_HPP
