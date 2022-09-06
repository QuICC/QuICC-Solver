/**
 * @file ISpectralKernel.hpp
 * @brief Base building block for the implementation of a spectral kernel 
 */

#ifndef QUICC_SPECTRAL_KERNEL_ISPECTRALKERNEL_HPP
#define QUICC_SPECTRAL_KERNEL_ISPECTRALKERNEL_HPP

// First include
//

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   /**
    * @brief Base building block for the implementation of a spectral kernel
    */
   class ISpectralKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit ISpectralKernel(const bool isComplex);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ISpectralKernel();

         /**
          * @brief Set the smart pointer to the scalar field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the scalar field
          */
         virtual void setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField);

         /**
          * @brief Set the smart pointer to the vector field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         virtual void setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField);

         /**
          * @brief Set resolution
          */
         void setResolution(SharedCResolution spRes);

         /**
          * @brief Compute the spectral kernel
          *
          * @param id   Component ID
          * @param i    Fastest index
          * @param j    Second index
          * @param k    Slowest index
          */
         virtual MHDVariant compute(const int i, const int j, const int k) const = 0;

         /**
          * @brief Apply kernel
          *
          * @param id   Component ID
          */
         virtual void apply();
         
      protected:
         /**
          * @brief Get scalar variable
          *
          * @param name Physical name of the field
          */
         const Framework::Selector::VariantSharedScalarVariable& scalar(std::size_t name) const;

         /**
          * @brief Get vector variable
          *
          * @param name Physical name of the field
          */
         const Framework::Selector::VariantSharedVectorVariable& vector(std::size_t name) const;

         /**
          * @brief Get shared resolution
          */
         SharedCResolution spRes() const;

         /**
          * @brief Get shared resolution
          */
         const Resolution& res() const;

         /**
          * @brief Spectral values are complex?
          */
         bool mIsComplex;

         /**
          * @brief Map of name and pointer for the scalar variables
          */
         std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>  mScalars;

         /**
          * @brief Map of name and pointer for the vector variables
          */
         std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>  mVectors;

      private:
         /**
          * @brief Shared resolution
          */
         SharedCResolution mspRes;

   };

   /// Typedef for a smart ISpectralKernel
   typedef std::shared_ptr<ISpectralKernel> SharedISpectralKernel;
   
}
}
}

#endif // QUICC_SPECTRAL_KERNEL_ISPECTRALKERNEL_HPP
