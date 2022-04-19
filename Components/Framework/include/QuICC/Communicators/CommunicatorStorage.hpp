/**
 * @file CommunicatorStorage.hpp
 * @brief Implementation of the storage providers for the communicators
 */

#ifndef QUICC_PARALLEL_COMMUNICATORSTORAGE_HPP
#define QUICC_PARALLEL_COMMUNICATORSTORAGE_HPP

// Configuration includes
//

// Configuration includes
//

// System includes
//
#include <vector>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/StorageProviders/DynamicPairProvider.hpp"
#include "QuICC/Communicators/Converters/IConverter.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Specialisation of the implementation of the storage providers for the communicators
    */
   class CommunicatorStorage
   {
      public:
         /// Typedef for storage type
         typedef DynamicPairProvider StorageType;

         /// Typedef for Shared setup type
         typedef typename StorageType::SharedFwdSetupType SharedFwdSetupType;

         /// Typedef for Shared setup type
         typedef typename StorageType::SharedBwdSetupType SharedBwdSetupType;

         /**
         * @brief Constructor
         */
         CommunicatorStorage();

         /**
         * @brief Destructor
         */
         virtual ~CommunicatorStorage();

         /**
          * @brief Check if converter exits
          */
         bool hasConverter(Dimensions::Transform::Id id) const;

         /**
          * @brief Add storage
          */
         void addStorage(Dimensions::Transform::Id id);

         /**
          * @brief Add converter
          */
         void addConverter(Dimensions::Transform::Id id, std::shared_ptr<IConverter> spConv);

         /**
          * @brief Get/Set storage provider
          */
         StorageType& storage(const Dimensions::Transform::Id id);

         /**
          * @brief Get/Set data converter
          */
         IConverter& converter(const Dimensions::Transform::Id id);

         /**
          * @brief Get/Set storage provider
          */
         template <Dimensions::Transform::Id TID> StorageType& storage();

         /**
          * @brief Get/Set data converter
          */
         template <Dimensions::Transform::Id TID> IConverter& converter();

         /**
          * @brief Get/Set last storage provider
          */
         StorageType& lastStorage();

      protected:

      private:
         /**
          * @brief Last storage ID
          */
         Dimensions::Transform::Id mLastStorageId;

         /**
          * @brief Storage providers
          */
         std::map<Dimensions::Transform::Id,StorageType>  mStorages;

         /**
          * @brief Storage providers
          */
         std::map<Dimensions::Transform::Id,std::shared_ptr<IConverter> >  mConverters;
   };

   inline CommunicatorStorage::StorageType& CommunicatorStorage::storage(const Dimensions::Transform::Id id)
   {
      return this->mStorages.find(id)->second;
   }

   inline IConverter& CommunicatorStorage::converter(const Dimensions::Transform::Id id)
   {
      return *this->mConverters.find(id)->second;
   }

   template <Dimensions::Transform::Id TID>
      typename CommunicatorStorage::StorageType& CommunicatorStorage::storage()
   {
      return this->mStorages.find(TID)->second;
   }

   template <Dimensions::Transform::Id TID>
      IConverter& CommunicatorStorage::converter()
   {
      return *this->mConverters.find(TID)->second;
   }

}
}

#endif // QUICC_PARALLEL_COMMUNICATORSTORAGE_HPP
