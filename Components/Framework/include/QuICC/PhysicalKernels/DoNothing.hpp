/**
 * @file DoNothing.hpp
 * @brief Trivial physical kernel that does nothing
 */

#ifndef QUICC_PHYSICAL_KERNEL_DONOTHING_HPP
#define QUICC_PHYSICAL_KERNEL_DONOTHING_HPP

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

   /**
    * @brief Trivial kernel that does nothing
    */
   class DoNothing: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit DoNothing();

         /**
          * @brief Simple empty destructor
          */
         virtual ~DoNothing();

         /**
          * @brief Set the smart pointer to the scalar field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the scalar field
          */
         virtual void setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField) override;

         /**
          * @brief Set the smart pointer to the vector field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         virtual void setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField) override;

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

   };

   /// Typedef for a smart DoNothing
   typedef std::shared_ptr<DoNothing> SharedDoNothing;
   
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_DONOTHING_HPP
