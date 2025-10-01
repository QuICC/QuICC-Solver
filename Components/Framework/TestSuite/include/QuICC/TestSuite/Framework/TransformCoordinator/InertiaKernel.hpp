/**
 * @file InertiaKernel.hpp
 * @brief Physical kernel for the Momentum nonlinear kernel
 */

#ifndef QUICC_PHYSICAL_INERTIAKERNEL_HPP
#define QUICC_PHYSICAL_INERTIAKERNEL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   /**
    * @brief Physical kernel for the Momentum nonlinear kernel
    */
   class InertiaKernel: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit InertiaKernel();

         /**
          * @brief Simple empty destructor
          */
         virtual ~InertiaKernel() = default;

         /**
          * @brief Set the physical mesh on which kernel is working
          */
         virtual void setMesh(std::shared_ptr<std::vector<Array> > spMesh) override;

         /**
          * @brief Set the smart pointer to the velocitt field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         virtual void setVelocity(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField);

         /**
          * @brief Initialize kernel
          */
         void init(const MHDFloat inertia);

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
          * @brief Name ID of the unknown
          */
         std::size_t mName;

         /**
          * @brief Name ID of the temperature field
          */
         std::size_t mTempName;

         /**
          * @brief Scaling constant for inertial term
          */
         MHDFloat mInertia;

         /**
          * @brief Storage for the radial grid values (if required)
          */
         Array mRadius;

   };

}
}
}

#endif // QUICC_PHYSICAL_INERTIAKERNEL_HPP
