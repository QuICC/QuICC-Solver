/**
 * @file DynamicRZProvider.hpp
 * @brief Implementation of a combined read and complex data storage provider adapting its size dynamically.
 */

#ifndef QUICC_DYNAMICRZPROVIDER_HPP
#define QUICC_DYNAMICRZPROVIDER_HPP

// System includes
//
#include <cassert>
#include <stdexcept>
#include <variant>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/StorageProviders/DynamicStorageProvider.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a combined real and complex data storage provider adapting its size dynamically.
    */
   class DynamicRZProvider
   {
      public:
         /// Typedef for real data storage provider
         typedef DynamicStorageProvider<MHDFloat>  RealStorage;

         /// Typedef for complex data storage provider
         typedef DynamicStorageProvider<MHDComplex>  ComplexStorage;

         /// Typedef for variant of both storage types
         typedef std::variant<RealStorage::DataType *, ComplexStorage::DataType *> VariantDataPointer;

         /// Typedef for setup
         typedef RealStorage::SharedSetupType SharedSetupType;

         /**
          * @brief Constructor
          */
         DynamicRZProvider();

         /**
          * @brief Destructor
          */
         ~DynamicRZProvider();

         /**
          * @brief Initialise the provider
          *
          * @param spSetup Setup for datatype
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Provide real or complex temporary storage
          */
         void provide(VariantDataPointer& pTmp);

         /**
          * @brief Provide real temporary storage
          */
         void provide(RealStorage::DataType *& pTmp);

         /**
          * @brief Provide complex temporary storage
          */
         void provide(ComplexStorage::DataType *& pTmp);

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void hold(VariantDataPointer& pTmp);

         /**
          * @brief Hold real storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void hold(RealStorage::DataType &tmp);

         /**
          * @brief Hold complex storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void hold(ComplexStorage::DataType &tmp);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used map)
          * @param id  ID used as key in map
          */
         void hold(VariantDataPointer& pTmp, const int id);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used map)
          * @param id  ID used as key in map
          */
         void hold(RealStorage::DataType &tmp, const int id);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to free (put in used map)
          * @param id  ID used as key in map
          */
         void hold(ComplexStorage::DataType &tmp, const int id);

         /**
          * @brief Recover the hold temporary storage used previously
          */
         void recover(VariantDataPointer& pTmp);

         /**
          * @brief Recover the hold real temporary storage used previously
          */
         void recover(RealStorage::DataType *& pTmp);

         /**
          * @brief Recover the hold complex temporary storage used previously
          */
         void recover(ComplexStorage::DataType *& pTmp);

         /**
          * @brief Recover the hold temporary storage with specific ID used previously
          */
         void recover(VariantDataPointer& pTmp, const int id);

         /**
          * @brief Recover the hold real temporary storage with specific ID used previously
          */
         void recover(RealStorage::DataType *& pTmp, const int id);

         /**
          * @brief Recover the hold complex temporary storage with specifig ID used previously
          */
         void recover(ComplexStorage::DataType *& pTmp, const int id);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void free(VariantDataPointer& pTmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void free(RealStorage::DataType &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void free(ComplexStorage::DataType &tmp);

         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Real data storage pool
          */
         RealStorage mReal;

         /**
          * @brief Complex data storage pool
          */
         ComplexStorage mComplex;
   };

   namespace internal {
      /// Visitor for free call
      struct FreeVisitor {
         FreeVisitor(DynamicRZProvider& p): storage(p) {};
         template <typename T> inline void operator()(T *& arg){storage.free(*arg);};
         DynamicRZProvider& storage;
      };
      /// Visitor for hold call
      struct HoldVisitor {
         HoldVisitor(DynamicRZProvider& p): storage(p) {};
         template <typename T> inline void operator()(T *& arg){storage.hold(*arg);};
         DynamicRZProvider& storage;
      };
      /// Visitor for hold call
      struct HoldIdVisitor {
         HoldIdVisitor(DynamicRZProvider& p, const int id): storage(p), id(id) {};
         template <typename T> inline void operator()(T *& arg){storage.hold(*arg, id);};
         DynamicRZProvider& storage;
         const int id;
      };
      /// Visitor for recover call
      struct RecoverVisitor {
         RecoverVisitor(DynamicRZProvider& p): storage(p) {};
         template <typename T> inline void operator()(T *& arg){storage.recover(arg);};
         DynamicRZProvider& storage;
      };
      /// Visitor for recover ID call
      struct RecoverIdVisitor {
         RecoverIdVisitor(DynamicRZProvider& p, const int id): storage(p), id(id) {};
         template <typename T> inline void operator()(T *& arg){storage.recover(arg, id);};
         DynamicRZProvider& storage;
         const int id;
      };
      /// Visitor for provide call
      struct ProvideVisitor {
         ProvideVisitor(DynamicRZProvider& p): storage(p) {};
         template <typename T> void operator()(T *& arg){storage.provide(arg);};
         DynamicRZProvider& storage;
      };
   }

   inline void DynamicRZProvider::free(VariantDataPointer &pTmp)
   {
      std::visit(internal::FreeVisitor(*this), pTmp);
   }

   inline void DynamicRZProvider::free(RealStorage::DataType &tmp)
   {
      this->mReal.free(tmp);
   }

   inline void DynamicRZProvider::free(ComplexStorage::DataType &tmp)
   {
      this->mComplex.free(tmp);
   }

   inline void DynamicRZProvider::hold(VariantDataPointer &pTmp)
   {
      std::visit(internal::HoldVisitor(*this), pTmp);
   }

   inline void DynamicRZProvider::hold(RealStorage::DataType &tmp)
   {
      this->mReal.hold(tmp);
   }

   inline void DynamicRZProvider::hold(ComplexStorage::DataType &tmp)
   {
      this->mComplex.hold(tmp);
   }

   inline void DynamicRZProvider::hold(VariantDataPointer &pTmp, const int id)
   {
      std::visit(internal::HoldIdVisitor(*this,id), pTmp);
   }

   inline void DynamicRZProvider::hold(RealStorage::DataType &tmp, const int id)
   {
      this->mReal.hold(tmp,id);
   }

   inline void DynamicRZProvider::hold(ComplexStorage::DataType &tmp, const int id)
   {
      this->mComplex.hold(tmp,id);
   }

   inline void DynamicRZProvider::provide(VariantDataPointer& pTmp)
   {
      std::visit(internal::ProvideVisitor(*this), pTmp);
   }

   inline void DynamicRZProvider::provide(RealStorage::DataType *& pTmp)
   {
      this->mReal.provide(pTmp);
   }

   inline void DynamicRZProvider::provide(ComplexStorage::DataType *& pTmp)
   {
      this->mComplex.provide(pTmp);
   }

   inline void DynamicRZProvider::recover(VariantDataPointer& pTmp)
   {
      std::visit(internal::RecoverVisitor(*this), pTmp);
   }

   inline void DynamicRZProvider::recover(RealStorage::DataType *& pTmp)
   {
      this->mReal.recover(pTmp);
   }

   inline void DynamicRZProvider::recover(ComplexStorage::DataType *& pTmp)
   {
      this->mComplex.recover(pTmp);
   }

   inline void DynamicRZProvider::recover(VariantDataPointer& pTmp ,const int id)
   {
      std::visit(internal::RecoverIdVisitor(*this,id), pTmp);
   }

   inline void DynamicRZProvider::recover(RealStorage::DataType *& pTmp ,const int id)
   {
      this->mReal.recover(pTmp,id);
   }

   inline void DynamicRZProvider::recover(ComplexStorage::DataType *& pTmp, const int id)
   {
      this->mComplex.recover(pTmp,id);
   }

}

#endif // QUICC_DYNAMICRZPROVIDER_HPP
