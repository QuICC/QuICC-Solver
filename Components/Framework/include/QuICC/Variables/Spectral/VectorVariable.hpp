/**
 * @file VectorVariable.hpp
 * @brief Implementation of vector variable
 */

#ifndef QUICC_DATATYPES_VECTORVARIABLE_HPP
#define QUICC_DATATYPES_VECTORVARIABLE_HPP

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
#include "QuICC/Variables/Physical/VectorPhysicalVariable.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation of vector variable
    */
   template <typename TSScalar, typename TPScalar> class VectorVariable: public VectorPhysicalVariable<TPScalar>
   {
      public:
         /**
          * @brief Constructs the underlying physical and spectral fields
          *
          * @param spRes Resolution information
          */
         VectorVariable(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~VectorVariable();

         /**
          * @brief Get spectral vector field (perturbation part)
          */
         const VectorField<TSScalar,FieldComponents::Spectral::Id>&  perturbation() const;

         /**
          * @brief Set spectral vector field (perturbation part)
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>&  rPerturbation();

         /**
          * @brief Get spectral vector field (total field)
          */
         const VectorField<TSScalar,FieldComponents::Spectral::Id>&  total() const;

         /**
          * @brief Set spectral vector field (total field)
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>&  rTotal();

         /**
          * @brief initialise to zeros
          */
         void setZeros();

         /**
          * @brief Initialise the spectral values storage
          */
         void initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps);

         /**
         * @brief Get the memory requirements
         */
         virtual MHDFloat requiredStorage() const;

         /**
         * @brief Profile the memory requirements
         */
         virtual void profileStorage() const;

      protected:

         /**
          * @brief Spectral vector of the field
          */
         std::shared_ptr<VectorField<TSScalar,FieldComponents::Spectral::Id> > mspPerturbation;

      private:
   };

   template <typename TSScalar, typename TPScalar> inline const VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::perturbation() const
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::rPerturbation()
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline const VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::total() const
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::rTotal()
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> VectorVariable<TSScalar,TPScalar>::VectorVariable(SharedResolution spRes)
      : VectorPhysicalVariable<TPScalar>(spRes)
   {
   }

   template <typename TSScalar, typename TPScalar> VectorVariable<TSScalar,TPScalar>::~VectorVariable()
   {
   }

   template <typename TSScalar, typename TPScalar> void VectorVariable<TSScalar,TPScalar>::setZeros()
   {
      // initialise the physical components to zero
      VectorPhysicalVariable<TPScalar>::setZeros();

      // initialise vector field to zero
      this->rPerturbation().setZeros();
   }

   template <typename TSScalar, typename TPScalar> void VectorVariable<TSScalar,TPScalar>::initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps)
   {
      // Safety assert
      assert(! this->mspPerturbation);

      std::map<FieldComponents::Spectral::Id,bool> map;
      for(unsigned int i = 0; i < comps.size(); i++)
      {
         map.insert(std::make_pair(comps.at(i), true));
      }

      this->mspPerturbation = std::make_shared<VectorField<TSScalar,FieldComponents::Spectral::Id> >(this->res().spSpectralSetup(), map);
   }

   template <typename TSScalar, typename TPScalar> MHDFloat VectorVariable<TSScalar,TPScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += VectorPhysicalVariable<TPScalar>::requiredStorage();

      mem += this->mspPerturbation->requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   template <typename TSScalar, typename TPScalar> void VectorVariable<TSScalar,TPScalar>::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = VectorPhysicalVariable<TPScalar>::requiredStorage();
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLESPHYS, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLES, mem);

      mem = this->mspPerturbation->requiredStorage();
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLESSPEC, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLES, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}

#endif // QUICC_DATATYPES_VECTORVARIABLE_HPP
