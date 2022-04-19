/**
 * @file MakeConstant.hpp
 * @brief Trivial kernel to make field constant
 */

#ifndef QUICC_PHYSICAL_KERNEL_MAKECONSTANT_HPP
#define QUICC_PHYSICAL_KERNEL_MAKECONSTANT_HPP

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
   class MakeConstant: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit MakeConstant();

         /**
          * @brief Simple empty destructor
          */
         virtual ~MakeConstant();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDFloat value);

         /**
          * @brief Compute the physical kernel
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const;
         
      protected:

      private:
         /**
          * @brief Value to set the field to
          */
         MHDFloat mValue;

   };

   /// Typedef for a smart MakeConstant
   typedef std::shared_ptr<MakeConstant> SharedMakeConstant;
   
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_MAKECONSTANT_HPP
