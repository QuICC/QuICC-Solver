/**
 * @file IPhysicalKernel.hpp
 * @brief Base building block for the implementation of a physical kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_IPHYSICALKERNEL_HPP
#define QUICC_PHYSICAL_KERNEL_IPHYSICALKERNEL_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   /**
    * @brief Base building block for the implementation of a physical kernel
    */
   class IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit IPhysicalKernel();

         /**
          * @brief Simple empty destructor
          */
         virtual ~IPhysicalKernel();

         /**
          * @brief Set the physical mesh on which kernel is working
          */
         virtual void setMesh(std::shared_ptr<std::vector<Array> > spMesh);

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
          * @brief Compute the physical kernel
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void compute(Framework::Selector::RealScalarField& rNLComp, FieldComponents::Physical::Id id) const = 0;

         /**
          * @brief Compute the physical kernel (required for variant, probably never implemented)
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void compute(Framework::Selector::ComplexScalarField& rNLComp, FieldComponents::Physical::Id id) const;

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
          * @brief Mesh used by the kernel
          */
         std::shared_ptr<std::vector<Array> > mspMesh;

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

   /// Typedef for a smart IPhysicalKernel
   typedef std::shared_ptr<IPhysicalKernel> SharedIPhysicalKernel;

}
}
}

#endif // QUICC_PHYSICAL_KERNEL_IPHYSICALKERNEL_HPP
