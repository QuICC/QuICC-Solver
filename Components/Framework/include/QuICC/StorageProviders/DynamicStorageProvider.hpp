/** 
 * @file DynamicStorageProvider.hpp
 * @brief Templated implementation of a data storage provider adapting its size dynamically.
 */

#ifndef QUICC_DYNAMICSTORAGEPROVIDER_HPP
#define QUICC_DYNAMICSTORAGEPROVIDER_HPP

// System includes
//
#include <cassert>
#include <queue>
#include <list>
#include <stdexcept>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

   /**
    * @brief Templated implementation of a data storage provider adapting its size dynamically.
    *
    * \tparam T Data type
    */
   template <typename T> class DynamicStorageProvider
   {
      public:
         /// Typedef for storage data type
         typedef Framework::Selector::ScalarField<T>  DataType;

         /// Typedef for storage data type
         typedef typename DataType::SharedSetupType SharedSetupType;

         /**
          * @brief Constructor
          */
         DynamicStorageProvider();

         /**
          * @brief Destructor
          */
         ~DynamicStorageProvider();

         /**
          * @brief Initialise the provider with one unit of each
          *
          * \mhdBug Should ultimatively be removed or at least simplified
          *
          * @param spSetup Setup for storage
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Provide temporary storage
          */
         void provide(DataType *& pTmp);

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void hold(DataType &tmp);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used map)
          * @param id  ID used as key in map
          */
         void hold(DataType &tmp, const int id);

         /**
          * @brief Recover the unfreed temporary storage used previously
          */
         void recover(DataType *& pTmp);

         /**
          * @brief Recover the unfreed temporary storage with specific ID used previously
          */
         void recover(DataType *& pTmp, const int id);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void free(DataType &tmp);

         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
         
      protected:

      private:
         /**
          * @brief Reset all the pointer queues
          */
         void reset();

         /**
          * @brief Add storage
          */
         void addStorage();

         /**
          * @brief Maximum size of dynamic pool
          */
         const typename std::list<DataType>::size_type MAX_SIZE;

         /**
          * @brief Shared setup
          */
         SharedSetupType mspSetup;

         /**
          * @brief Available temporary storage queue
          */
         std::queue<DataType *> mAvailable;

         /**
          * @brief In-use temporary storage queue
          */
         std::queue<DataType *> mInUse;

         /**
          * @brief In-use temporary storage map
          */
         std::map<int,DataType *> mInUseMap;

         /**
          * @brief Storage pool
          */
         std::list<DataType>  mPool;
   };

   template <typename T> inline void DynamicStorageProvider<T>::free(DataType &tmp)
   {
      // Check if data is part of pool
      assert(std::any_of(this->mPool.begin(), this->mPool.end(),
            [&](const DataType& d)
            {
               return &d == &tmp;
            }));

      this->mAvailable.push(&tmp);
   }

   template <typename T> inline void DynamicStorageProvider<T>::hold(DataType &tmp)
   {
      this->mInUse.push(&tmp);
   }

   template <typename T> inline void DynamicStorageProvider<T>::hold(DataType &tmp, const int id)
   {
      this->mInUseMap.insert(std::make_pair(id, &tmp));
   }

   template <typename T> void DynamicStorageProvider<T>::provide(DataType *& pTmp)
   {
      // Check if storage is avaiable
      if(this->mAvailable.empty())
      {
         this->addStorage();
      }

      // Assert for non empty storage
      assert(!this->mAvailable.empty());

      // Get pointer from available queue
      pTmp = this->mAvailable.front();

      // Remove pointer from availabe queue
      this->mAvailable.pop();
   }

   template <typename T> void DynamicStorageProvider<T>::recover(DataType *& pTmp)
   {
      // Add in an assert for non empty queue
      assert(this->mInUse.size() > 0);

      // Get pointer from used queue
      pTmp = this->mInUse.front();

      // Remove pointer from used queue
      this->mInUse.pop();
   }

   template <typename T> void DynamicStorageProvider<T>::recover(DataType *& pTmp, const int id)
   {
      // Add in an assert for non empty map
      assert(this->mInUseMap.count(id) > 0);

      // Get pointer from used map
      pTmp = this->mInUseMap.find(id)->second;

      // Remove pointer from used map
      this->mInUseMap.erase(id);
   }

   template <typename T> DynamicStorageProvider<T>::DynamicStorageProvider()
      : MAX_SIZE(20)
   {
   }

   template <typename T> DynamicStorageProvider<T>::~DynamicStorageProvider()
   {
   }

   template <typename T> void DynamicStorageProvider<T>::init(SharedSetupType spSetup)
   {
      // Set the setup
      this->mspSetup = spSetup;

      // Reset all queues and pool
      this->reset();

      // Add first storage
      this->addStorage();
   }

   template <typename T> void DynamicStorageProvider<T>::addStorage()
   {
      // Avoid growing to big (error in usage most likely)
      if(this->mPool.size() > MAX_SIZE)
      {
         throw std::logic_error("Dynamic storage provider exceeded MAX_SIZE!");
      }

      // Create storage scalar
      this->mPool.push_back(DataType(this->mspSetup));

      // Add storage to queue
      this->mAvailable.push(&(this->mPool.back()));
   }

   template <typename T> void DynamicStorageProvider<T>::reset()
   {
      // Make sure available queue is empty
      if(!this->mAvailable.empty())
      {
         while(!this->mAvailable.empty()) this->mAvailable.pop();
      }

      // Make sure in-use queue is empty
      if(!this->mInUse.empty())
      {
         while(!this->mInUse.empty()) this->mInUse.pop();
      }

      // Make sure in-use map is empty
      this->mInUseMap.clear();

      // Empty the data pool
      this->mPool.clear();
   }

   template <typename T> MHDFloat  DynamicStorageProvider<T>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

   #ifdef QUICC_STORAGEPROFILE
      if(this->mPool.size() > 0)
      {
         mem += this->mPool.size()*this->mPool.front().requiredStorage();
      }
   #endif // QUICC_STORAGEPROFILE

      return mem;
   }

}

#endif // QUICC_DYNAMICSTORAGEPROVIDER_HPP
