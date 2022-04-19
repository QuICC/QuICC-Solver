/**
 * @file Passthrough.hpp
 * @brief Trivial passthrough kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_PASSTHROUGH_HPP
#define QUICC_PHYSICAL_KERNEL_PASSTHROUGH_HPP

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
    * @brief Trivial passthrough kernel
    */
   class Passthrough: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit Passthrough();

         /**
          * @brief Simple empty destructor
          */
         virtual ~Passthrough();

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
          * @brief Compute the physical kernel
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const;
         
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

   /// Typedef for a smart Passthrough
   typedef std::shared_ptr<Passthrough> SharedPassthrough;
   
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_PASSTHROUGH_HPP
