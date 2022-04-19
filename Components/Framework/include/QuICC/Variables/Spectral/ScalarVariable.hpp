/**
 * @file ScalarVariable.hpp
 * @brief Implementation of scalar field variable
 */

#ifndef QUICC_DATATYPES_SCALARVARIABLE_HPP
#define QUICC_DATATYPES_SCALARVARIABLE_HPP

// Configuration includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Variables/Physical/ScalarPhysicalVariable.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation of scalar field variable
    */
   template <typename TSScalar, typename TPScalar> class ScalarVariable: public ScalarPhysicalVariable<TPScalar>
   {
      public:
         /**
          * @brief Constructs the underlying physical and spectral fields
          *
          * @param spRes Resolution information
          */
         ScalarVariable(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~ScalarVariable();

         /**
          * @brief Get scalar field (perturbation part)
          */
         const TSScalar&  perturbation() const;

         /**
          * @brief Set scalar field (perturbation part)
          */
         TSScalar&  rPerturbation();

         /**
          * @brief Get scalar field (total field)
          */
         const TSScalar&  total() const;

         /**
          * @brief Set scalar field (total field)
          */
         TSScalar&  rTotal();

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
          * @brief Spectral scalar field
          */
         TSScalar    mPerturbation;

      private:
   };

   template <typename TSScalar, typename TPScalar> inline const TSScalar& ScalarVariable<TSScalar,TPScalar>::perturbation() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline TSScalar& ScalarVariable<TSScalar,TPScalar>::rPerturbation()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline const TSScalar& ScalarVariable<TSScalar,TPScalar>::total() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline TSScalar& ScalarVariable<TSScalar,TPScalar>::rTotal()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> ScalarVariable<TSScalar,TPScalar>::ScalarVariable(SharedResolution spRes)
      : ScalarPhysicalVariable<TPScalar>(spRes), mPerturbation(spRes->spSpectralSetup())
   {
   }

   template <typename TSScalar, typename TPScalar> ScalarVariable<TSScalar,TPScalar>::~ScalarVariable()
   {
   }

   template <typename TSScalar, typename TPScalar> void ScalarVariable<TSScalar,TPScalar>::setZeros()
   {
      // initialise the physical components to zero
      ScalarPhysicalVariable<TPScalar>::setZeros();

      // initialise scalar field to zero
      this->rPerturbation().setZeros();
   }

   template <typename TSScalar, typename TPScalar> void ScalarVariable<TSScalar,TPScalar>::initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps)
   {
   }

   template <typename TSScalar, typename TPScalar> MHDFloat ScalarVariable<TSScalar,TPScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += ScalarPhysicalVariable<TPScalar>::requiredStorage();

      mem += this->mPerturbation.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   template <typename TSScalar, typename TPScalar> void ScalarVariable<TSScalar,TPScalar>::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = ScalarPhysicalVariable<TPScalar>::requiredStorage();
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLESPHYS, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLES, mem);

      mem = this->mPerturbation.requiredStorage();
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLESSPEC, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::VARIABLES, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}

#endif // QUICC_DATATYPES_SCALARVARIABLE_HPP
