/**
 * @file VerticalField.hpp
 * @brief Physical kernel for the vertical component of spherical field
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERICAL_VERTICALFIELD_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERICAL_VERTICALFIELD_HPP

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
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Spherical {

   /**
    * @brief Physical kernel for the vertical component of spherical field
    */
   class VerticalField: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit VerticalField();

         /**
          * @brief Simple empty destructor
          */
         virtual ~VerticalField();

         /**
          * @brief Set the physical mesh on which kernel is working
          */
         virtual void setMesh(std::shared_ptr<std::vector<Array> > spMesh) override;

         /**
          * @brief Set the smart pointer to the vector field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         void setVector(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField);

         /**
          * @brief Initialize kernel
          */
         void init(FieldType::Id type, const MHDFloat transport);

         /**
          * @brief Compute the physical kernel
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const override;
         
      protected:
         /**
          * @brief Get name ID of the unknown
          */
         std::size_t name() const;

      private:
         /**
          * @brief Name ID
          */
         std::size_t mName;

         /**
          * @brief Storage for output field flag
          */
         FieldType::Id mFieldType;

         /**
          * @brief Scaling constant
          */
         MHDFloat mScale;

         /**
          * @brief Storage for the cos(theta) grid values (if required)
          */
         Array mCosTheta;

         /**
          * @brief Storage for the sin(theta) grid values (if required)
          */
         Array mSinTheta;

   };

   /// Typedef for a smart VerticalField
   typedef std::shared_ptr<VerticalField> SharedVerticalField;
   
}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERICAL_VERTICALFIELD_HPP
