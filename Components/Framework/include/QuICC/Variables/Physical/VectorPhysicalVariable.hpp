/**
 * @file VectorPhysicalVariable.hpp
 * @brief Base of the implementation of the physical components of a vector variable
 */

#ifndef QUICC_DATATYPES_VECTORPHYSICALVARIABLE_HPP
#define QUICC_DATATYPES_VECTORPHYSICALVARIABLE_HPP

// Configuration includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/TensorFields/SymmetricTensorField.hpp"
#include "QuICC/Variables/VariableBase.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Base of the implementation of the physical components of a vector variable
    */
   template <typename TScalar> class VectorPhysicalVariable : public VariableBase
   {
      public:
         /**
         * @brief Constructor
         *
         * There is only very little done directly. The different fields have to be
         * initialised after construction.
         *
         * @param spRes Resolution information
         */
         VectorPhysicalVariable(SharedResolution spRes);

         /**
         * @brief Simple empty destructor
         */
         virtual ~VectorPhysicalVariable();

         /**
          * @brief Get the physical field values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   phys() const;

         /**
          * @brief Set the physical field values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rPhys();

         /**
          * @brief Get the physical vector gradient values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   grad(FieldComponents::Spectral::Id id) const;

         /**
          * @brief Set the physical vector gradient values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rGrad(FieldComponents::Spectral::Id id);

         /**
          * @brief Get the physical curl values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   curl() const;

         /**
          * @brief Set the physical curl values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rCurl();

         /**
          * @brief Get the physical vector 2nd order gradient values
          */
         const SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&   grad2(FieldComponents::Spectral::Id id) const;

         /**
          * @brief Set the physical vector 2nd order gradient values
          */
         SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&   rGrad2(FieldComponents::Spectral::Id id);

         /**
          * @brief Initialise to zero
          */
         void setZeros();

         /**
          * @brief Initialise the physical values storage
          */
         void initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical gradient storage
          */
         void initPhysicalGradient(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical curl storage
          */
         void initPhysicalCurl(const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical 2nd order gradient storage
          */
         void initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& comps);

         /**
          * @brief Check if variable has physical data setup
          */
         bool hasPhys() const;

         /**
          * @brief Check if variable has vector gradient data setup
          */
         bool hasGrad() const;

         /**
          * @brief Check if variable has curl data setup
          */
         bool hasCurl() const;

         /**
          * @brief Check if variable has vector 2nd order gradient data setup
          */
         bool hasGrad2() const;

         /**
         * @brief Get the memory requirements
         */
         virtual MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Smart pointer for the physical space field values
          */
         std::shared_ptr<VectorField<TScalar,FieldComponents::Physical::Id> > mspPhys;

         /**
          * @brief Smart pointer for the physical vector gradient values
          */
         std::map<FieldComponents::Spectral::Id,std::shared_ptr<VectorField<TScalar,FieldComponents::Physical::Id> > > mVGrad;

         /**
          * @brief Smart pointer for the physical curl values
          */
         std::shared_ptr<VectorField<TScalar,FieldComponents::Physical::Id> > mspCurl;

         /**
          * @brief Smart pointer for the physical vector 2nd order gradient values
          */
         std::map<FieldComponents::Spectral::Id,std::shared_ptr<SymmetricTensorField<TScalar,FieldComponents::Physical::Id> > > mVGrad2;
   };

   template <typename TScalar> inline bool VectorPhysicalVariable<TScalar>::hasPhys() const
   {
      return static_cast<bool>(this->mspPhys);
   }

   template <typename TScalar> inline bool VectorPhysicalVariable<TScalar>::hasGrad() const
   {
      return this->mVGrad.size();
   }

   template <typename TScalar> inline bool VectorPhysicalVariable<TScalar>::hasCurl() const
   {
      return static_cast<bool>(this->mspCurl);
   }

   template <typename TScalar> inline bool VectorPhysicalVariable<TScalar>::hasGrad2() const
   {
      return this->mVGrad2.size();
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::phys() const
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::rPhys()
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::grad(FieldComponents::Spectral::Id id) const
   {
      // Safety assertion
      assert(this->mVGrad.count(id));

      return *(this->mVGrad.find(id)->second);
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::rGrad(FieldComponents::Spectral::Id id)
   {
      // Safety assertion
      assert(this->mVGrad.count(id));

      return *(this->mVGrad.find(id)->second);
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::curl() const
   {
      // Safety assertion
      assert(this->mspCurl);

      return *this->mspCurl;
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::rCurl()
   {
      // Safety assertion
      assert(this->mspCurl);

      return *this->mspCurl;
   }

   template <typename TScalar> inline const SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::grad2(FieldComponents::Spectral::Id id) const
   {
      // Safety assertion
      assert(this->mVGrad2.count(id));

      return *(this->mVGrad2.find(id)->second);
   }

   template <typename TScalar> inline SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::rGrad2(FieldComponents::Spectral::Id id)
   {
      // Safety assertion
      assert(this->mVGrad2.count(id));

      return *(this->mVGrad2.find(id)->second);
   }

   template <typename TScalar> VectorPhysicalVariable<TScalar>::VectorPhysicalVariable(SharedResolution spRes)
      : VariableBase(spRes)
   {
   }

   template <typename TScalar> VectorPhysicalVariable<TScalar>::~VectorPhysicalVariable()
   {
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::setZeros()
   {
      // Initialise physical values to zero if required
      if(this->hasPhys())
      {
         this->rPhys().setZeros();
      }

      // Initialise vector gradient values to zero if required
      if(this->hasGrad())
      {
         for(auto it = this->mVGrad.begin(); it != this->mVGrad.end(); ++it)
         {
            it->second->setZeros();
         }
      }

      // Initialise curl values to zero if required
      if(this->hasCurl())
      {
         this->rCurl().setZeros();
      }

      // Initialise vector gradient values to zero if required
      if(this->mVGrad2.size() > 0)
      {
         for(auto it = this->mVGrad2.begin(); it != this->mVGrad2.end(); ++it)
         {
            it->second->setZeros();
         }
      }
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspPhys);

      this->mspPhys = std::make_shared<VectorField<TScalar,FieldComponents::Physical::Id> >(this->res().spPhysicalSetup(), comps);
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::initPhysicalGradient(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(this->mVGrad.count(id) == 0);

      // Create shared pointer
      auto spGrad = std::make_shared<VectorField<TScalar,FieldComponents::Physical::Id> >(this->res().spPhysicalSetup(), comps);

      // Insert into map
      this->mVGrad.insert(std::make_pair(id, spGrad));
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::initPhysicalCurl(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspCurl);

      this->mspCurl = std::make_shared<VectorField<TScalar,FieldComponents::Physical::Id> >(this->res().spPhysicalSetup(), comps);
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& comps)
   {
      // Safety assert
      assert(this->mVGrad2.count(id) == 0);

      // Create shared pointer
      auto spGrad2 = std::make_shared<SymmetricTensorField<TScalar,FieldComponents::Physical::Id> >(this->res().spPhysicalSetup(), comps);

      // Insert into map
      this->mVGrad2.insert(std::make_pair(id, spGrad2));
   }

   template <typename TScalar> MHDFloat VectorPhysicalVariable<TScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += VariableBase::requiredStorage();

      // Physical storage
      if(this->hasPhys())
      {
         mem += this->phys().requiredStorage();
      }

      // Physical vector gradient storage
      if(this->hasGrad())
      {
         for(auto it = this->mVGrad.cbegin(); it != this->mVGrad.cend(); ++it)
         {
            mem += it->second->requiredStorage();
         }
      }

      // Physical curl storage
      if(this->hasCurl())
      {
         mem += this->curl().requiredStorage();
      }

      // Physical vector 2nd order gradient storage
      if(this->hasGrad2())
      {
         for(auto it = this->mVGrad2.cbegin(); it != this->mVGrad2.cend(); ++it)
         {
            mem += it->second->requiredStorage();
         }
      }
#endif // QUICC_STORAGEPROFILE

      return mem;
   }
}
}

#endif // QUICC_DATATYPES_VECTORPHYSICALVARIABLE_HPP
