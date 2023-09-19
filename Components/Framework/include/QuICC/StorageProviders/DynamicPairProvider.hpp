/**
 * @file DynamicPairProvider.hpp
 * @brief Implementation of a data pair storage provider adapting its size dynamically.
 */

#ifndef QUICC_DYNAMICPAIRPROVIDER_HPP
#define QUICC_DYNAMICPAIRPROVIDER_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/StorageProviders/DynamicRZProvider.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a data pair storage provider adapting its size dynamically.
    */
   class DynamicPairProvider
   {
      public:
         /// Typedef for forward storage provider
         typedef DynamicRZProvider  FwdProvider;

         /// Typedef for forward storage provider
         typedef DynamicRZProvider  BwdProvider;

         /// Typedef for forward storage provider
         typedef FwdProvider::SharedSetupType SharedFwdSetupType;

         /// Typedef for forward storage provider
         typedef BwdProvider::SharedSetupType SharedBwdSetupType;

         /**
          * @brief Constructor
          */
         DynamicPairProvider();

         /**
          * @brief Destructor
          */
         ~DynamicPairProvider();

         /**
          * @brief Initialise the provider with one unit of each
          *
          * @param spSetupFwd Setup for TForward
          * @param spSetupBwd Setup for TBackward
          */
         void init(SharedFwdSetupType spSetupFwd, SharedBwdSetupType spSetupBwd);

         /**
          * @brief Provide temporary storage for forward transform
          */
         void  provideFwd(FwdProvider::VariantDataPointer & pTmp);

         /**
          * @brief Provide temporary storage for backward transform
          */
         void  provideBwd(BwdProvider::VariantDataPointer & pTmp);

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void  holdFwd(FwdProvider::VariantDataPointer &tmp);

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void  holdBwd(BwdProvider::VariantDataPointer &tmp);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used map)
          * @param id  ID used as key in map
          */
         void  holdFwd(FwdProvider::VariantDataPointer &tmp, const int id);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to free (put in used map)
          * @param id  ID used as key in map
          */
         void  holdBwd(BwdProvider::VariantDataPointer &tmp, const int id);

         /**
          * @brief Recover the unfreed temporary forward transform storage used previously
          */
         void  recoverFwd(FwdProvider::VariantDataPointer & pTmp);

         /**
          * @brief Recover the unfreed temporary backward transform storage used previously
          */
         void  recoverBwd(BwdProvider::VariantDataPointer & pTmp);

         /**
          * @brief Recover the unfreed temporary forward transform storage with specific ID used previously
          */
         void recoverFwd(FwdProvider::VariantDataPointer & pTmp, const int id);

         /**
          * @brief Recover the unfreed temporary backward transform storage with specifig ID used previously
          */
         void recoverBwd(BwdProvider::VariantDataPointer & pTmp, const int id);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void  freeFwd(FwdProvider::VariantDataPointer &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void  freeBwd(BwdProvider::VariantDataPointer &tmp);

         /**
          * @brief Provide temporary storage for forward transform
          */
         template <typename T> void  provideFwd(T & pTmp);

         /**
          * @brief Provide temporary storage for backward transform
          */
         template <typename T> void  provideBwd(T & pTmp);

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         template <typename T> void  holdFwd(T &tmp);

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to free (put in used queue)
          */
         template <typename T> void  holdBwd(T &tmp);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used map)
          * @param id  ID used as key in map
          */
         template <typename T>void  holdFwd(T &tmp, const int id);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to free (put in used map)
          * @param id  ID used as key in map
          */
         template <typename T> void  holdBwd(T &tmp, const int id);

         /**
          * @brief Recover the unfreed temporary forward transform storage used previously
          */
         template <typename T> void  recoverFwd(T & pTmp);

         /**
          * @brief Recover the unfreed temporary backward transform storage used previously
          */
         template <typename T> void  recoverBwd(T & pTmp);

         /**
          * @brief Recover the unfreed temporary forward transform storage with specific ID used previously
          */
         template <typename T> void recoverFwd(T & pTmp, const int id);

         /**
          * @brief Recover the unfreed temporary backward transform storage with specifig ID used previously
          */
         template <typename T> void recoverBwd(T & pTmp, const int id);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         template <typename T> void  freeFwd(T &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         template <typename T> void  freeBwd(T &tmp);

         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Storage provider for forward transform
          */
         FwdProvider mFwd;

         /**
          * @brief Storage provider for backward transform
          */
         FwdProvider mBwd;
   };

   inline void DynamicPairProvider::freeFwd(FwdProvider::VariantDataPointer& pTmp)
   {
      this->mFwd.free(pTmp);
   }

   inline void DynamicPairProvider::freeBwd(BwdProvider::VariantDataPointer& pTmp)
   {
      this->mBwd.free(pTmp);
   }

   inline void DynamicPairProvider::holdFwd(FwdProvider::VariantDataPointer& pTmp)
   {
      this->mFwd.hold(pTmp);
   }

   inline void DynamicPairProvider::holdBwd(BwdProvider::VariantDataPointer& pTmp)
   {
      this->mBwd.hold(pTmp);
   }

   inline void DynamicPairProvider::holdFwd(FwdProvider::VariantDataPointer& pTmp, const int id)
   {
      this->mFwd.hold(pTmp,id);
   }

   inline void DynamicPairProvider::holdBwd(BwdProvider::VariantDataPointer& pTmp, const int id)
   {
      this->mBwd.hold(pTmp, id);
   }

   inline void  DynamicPairProvider::provideFwd(FwdProvider::VariantDataPointer& pTmp)
   {
      this->mFwd.provide(pTmp);
   }

   inline void DynamicPairProvider::provideBwd(BwdProvider::VariantDataPointer& pTmp)
   {
      this->mBwd.provide(pTmp);
   }

   inline void DynamicPairProvider::recoverFwd(FwdProvider::VariantDataPointer& pTmp)
   {
      this->mFwd.recover(pTmp);
   }

   inline void DynamicPairProvider::recoverBwd(BwdProvider::VariantDataPointer& pTmp)
   {
      this->mBwd.recover(pTmp);
   }

   inline void DynamicPairProvider::recoverFwd(FwdProvider::VariantDataPointer& pTmp, const int id)
   {
      this->mFwd.recover(pTmp, id);
   }

   inline void DynamicPairProvider::recoverBwd(BwdProvider::VariantDataPointer& pTmp, const int id)
   {
      this->mBwd.recover(pTmp, id);
   }

   template <typename T> inline void DynamicPairProvider::freeFwd(T &tmp)
   {
      this->mFwd.free(tmp);
   }

   template <typename T> inline void DynamicPairProvider::freeBwd(T &tmp)
   {
      this->mBwd.free(tmp);
   }

   template <typename T> inline void DynamicPairProvider::holdFwd(T &tmp)
   {
      this->mFwd.hold(tmp);
   }

   template <typename T> inline void DynamicPairProvider::holdBwd(T &tmp)
   {
      this->mBwd.hold(tmp);
   }

   template <typename T> inline void DynamicPairProvider::holdFwd(T &tmp, const int id)
   {
      this->mFwd.hold(tmp,id);
   }

   template <typename T> inline void DynamicPairProvider::holdBwd(T &tmp, const int id)
   {
      this->mBwd.hold(tmp, id);
   }

   template <typename T> inline void  DynamicPairProvider::provideFwd(T &pTmp)
   {
      this->mFwd.provide(pTmp);
   }

   template <typename T> inline void DynamicPairProvider::provideBwd(T &pTmp)
   {
      this->mBwd.provide(pTmp);
   }

   template <typename T> inline void DynamicPairProvider::recoverFwd(T &pTmp)
   {
      this->mFwd.recover(pTmp);
   }

   template <typename T> inline void DynamicPairProvider::recoverBwd(T &pTmp)
   {
      this->mBwd.recover(pTmp);
   }

   template <typename T> inline void DynamicPairProvider::recoverFwd(T &pTmp, const int id)
   {
      this->mFwd.recover(pTmp, id);
   }

   template <typename T> inline void DynamicPairProvider::recoverBwd(T &pTmp, const int id)
   {
      this->mBwd.recover(pTmp, id);
   }

}

#endif // QUICC_DYNAMICPAIRPROVIDER_HPP
