/**
 * @file CommunicatorStorage.cpp
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/CommunicatorStorage.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   CommunicatorStorage::CommunicatorStorage()
      : mLastStorageId(Dimensions::Transform::TRA1D)
   {
   }

   CommunicatorStorage::~CommunicatorStorage()
   {
   }

   CommunicatorStorage::StorageType& CommunicatorStorage::lastStorage()
   {
      return this->mStorages.find(this->mLastStorageId)->second;
   }

   bool CommunicatorStorage::hasConverter(Dimensions::Transform::Id id) const
   {
      return (this->mConverters.count(id) == 1);
   }

   void CommunicatorStorage::addStorage(Dimensions::Transform::Id id)
   {
      this->mStorages.insert(std::make_pair(id,StorageType()));

      if(static_cast<int>(id) > static_cast<int>(this->mLastStorageId))
      {
         this->mLastStorageId = id;
      }
   }

   void CommunicatorStorage::addConverter(Dimensions::Transform::Id id, std::shared_ptr<IConverter> spConv)
   {
      this->mConverters.insert(std::make_pair(id,spConv));
   }

}
}
