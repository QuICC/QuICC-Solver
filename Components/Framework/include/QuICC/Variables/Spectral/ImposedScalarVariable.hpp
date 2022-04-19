/**
 * @file ImposedScalarVariable.hpp
 * @brief Implementation of scalar field variable with imposed component
 */

#ifndef QUICC_DATATYPES_IMPOSEDSCALARVARIABLE_HPP
#define QUICC_DATATYPES_IMPOSEDSCALARVARIABLE_HPP

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
#include "QuICC/Variables/Spectral/ScalarVariable.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation of scalar field variable with imposed component
    */
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> class ImposedScalarVariable: public ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>
   {
      public:
         /**
          * @brief Constructs the underlying physical and spectral fields
          *
          * @param spRes Resolution information
          */
         ImposedScalarVariable(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~ImposedScalarVariable();

         /**
          * @brief Get scalar field (total field)
          */
         const TSScalar&  total() const;

         /**
          * @brief Set scalar field (total field)
          */
         TSScalar&  rTotal();

         /**
          * @brief Get the imposed scalar field
          */
         const TSScalar&  imposed() const;

         /**
          * @brief Set the imposed scalar field
          */
         TSScalar&  rImposed();

         /**
          * @brief initialise to zeros
          */
         void setZeros();

         /**
         * @brief Get the memory requirements
         */
         virtual MHDFloat requiredStorage() const;

      protected:

         /**
          * @brief Imposed spectral scalar field
          */
         TSScalar    mImposed;

         /**
          * @brief Total scalar field (imposed + perturbation)
          */
         TSScalar    mTotal;

      private:
   };

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::total() const
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rTotal()
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::imposed() const
   {
      return this->mImposed;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rImposed()
   {
      return this->mImposed;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ImposedScalarVariabl<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>e::ImposedScalarVariable(SharedResolution spRes)
      : ScalarVariable(spRes), mImposed(spRes->spSpectralSetup()), mTotal(spRes->spSpectralSetup())
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ImposedScalarVariabl<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>e::~ImposedScalarVariable()
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> void ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::setZeros()
   {
      // initialise the pertubation component to zero
      ScalarVariable::setZeros();

      // initialise imposed scalar field to zero
      this->rImposed().setZeros();

      // initialise total scalar field to zero
      this->mTotal.setZeros();
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> MHDFloat ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += ScalarVariable::requiredStorage();

      mem += this->mImposed.requiredStorage();

      mem += this->mTotal.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}

#endif // QUICC_DATATYPES_IMPOSEDSCALARVARIABLE_HPP
