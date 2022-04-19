/**
 * @file ScalarPhysicalVariable.hpp
 * @brief Base of the implementation of the physical components of a scalar variable
 */

#ifndef QUICC_DATATYPES_SCALARPHYSICALVARIABLE_HPP
#define QUICC_DATATYPES_SCALARPHYSICALVARIABLE_HPP

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
    * @brief Base of the implementation of the physical components of a scalar variable
    */
   template <typename TScalar> class ScalarPhysicalVariable : public VariableBase
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
         ScalarPhysicalVariable(SharedResolution spRes);

         /**
         * @brief Simple empty destructor
         */
         virtual ~ScalarPhysicalVariable();

         /**
          * @brief Get the physical field values
          */
         const TScalar&   phys() const;

         /**
          * @brief Set the physical field values
          */
         TScalar&   rPhys();

         /**
          * @brief Get the physical gradient values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   grad() const;

         /**
          * @brief Set the physical gradient values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rGrad();

         /**
          * @brief Get the physical 2nd order gradient values
          */
         const SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&   grad2() const;

         /**
          * @brief Set the physical 2nd order gradient values
          */
         SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&   rGrad2();

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
          * @brief Initialise the physical 2nd order gradient storage
          */
         void initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& comps);

         /**
          * @brief Check if variable has physical data setup
          */
         bool hasPhys() const;

         /**
          * @brief Check if variable has gradient data setup
          */
         bool hasGrad() const;

         /**
          * @brief Check if variable has 2nd order gradient data setup
          */
         bool hasGrad2() const;

         /**
         * @brief Get the memory requirements
         */
         virtual MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Smart pointer for the physical field values
          */
         std::shared_ptr<TScalar> mspPhys;

         /**
          * @brief Smart pointer for the physical gradient values
          */
         std::shared_ptr<VectorField<TScalar,FieldComponents::Physical::Id> > mspGrad;

         /**
          * @brief Smart pointer for the physical 2nd order gradient values
          */
         std::shared_ptr<SymmetricTensorField<TScalar,FieldComponents::Physical::Id> > mspGrad2;
   };

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasPhys() const
   {
      return static_cast<bool>(this->mspPhys);
   }

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasGrad() const
   {
      return static_cast<bool>(this->mspGrad);
   }

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasGrad2() const
   {
      return static_cast<bool>(this->mspGrad2);
   }

   template <typename TScalar> inline const TScalar&  ScalarPhysicalVariable<TScalar>::phys() const
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline TScalar&  ScalarPhysicalVariable<TScalar>::rPhys()
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::grad() const
   {
      // Safety assertion
      assert(this->mspGrad);

      return *this->mspGrad;
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::rGrad()
   {
      // Safety assertion
      assert(this->mspGrad);

      return *this->mspGrad;
   }

   template <typename TScalar> inline const SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::grad2() const
   {
      // Safety assertion
      assert(this->mspGrad2);

      return *this->mspGrad2;
   }

   template <typename TScalar> inline SymmetricTensorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::rGrad2()
   {
      // Safety assertion
      assert(this->mspGrad2);

      return *this->mspGrad2;
   }

   template <typename TScalar> ScalarPhysicalVariable<TScalar>::ScalarPhysicalVariable(SharedResolution spRes)
      : VariableBase(spRes)
   {
   }

   template <typename TScalar> ScalarPhysicalVariable<TScalar>::~ScalarPhysicalVariable()
   {
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::setZeros()
   {
      // Initialise physical values to zero if required
      if(this->mspPhys)
      {
         this->rPhys().setZeros();
      }

      // Initialise gradient values to zero if required
      if(this->mspGrad)
      {
         this->rGrad().setZeros();
      }

      // Initialise vector gradient values to zero if required
      if(this->mspGrad2)
      {
         this->rGrad2().setZeros();
      }
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspPhys);

      this->mspPhys = std::make_shared<TScalar>(this->res().spPhysicalSetup());
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysicalGradient(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspGrad);
      assert(id == FieldComponents::Spectral::SCALAR);

      this->mspGrad = std::make_shared<VectorField<TScalar,FieldComponents::Physical::Id> >(this->res().spPhysicalSetup(), comps);
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<std::pair<FieldComponents::Physical::Id, FieldComponents::Physical::Id>,bool>& comps)
   {
      // Safety assert
      assert(! this->mspGrad2);
      assert(id == FieldComponents::Spectral::SCALAR);

      this->mspGrad2 = std::make_shared<SymmetricTensorField<TScalar,FieldComponents::Physical::Id> >(this->res().spPhysicalSetup(), comps);
   }

   template <typename TScalar> MHDFloat ScalarPhysicalVariable<TScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += VariableBase::requiredStorage();

      // Physical storage
      if(this->hasPhys())
      {
         mem += this->phys().requiredStorage();
      }

      // Physical gradient storage
      if(this->hasGrad())
      {
         mem += this->grad().requiredStorage();
      }

      // Physical 2nd order gradient storage
      if(this->hasGrad2())
      {
         mem += this->grad2().requiredStorage();
      }
#endif // QUICC_STORAGEPROFILE

      return mem;
   }
}
}

#endif // QUICC_DATATYPES_SCALARPHYSICALVARIABLE_HPP
