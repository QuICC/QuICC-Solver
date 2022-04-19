/**
 * @file UniformAxialField.hpp
 * @brief Uniform axial field generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_UNIFORMAXIALFIELD_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_UNIFORMAXIALFIELD_HPP

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

namespace Sphere {

   /**
    * @brief Uniform axial field generator kernel
    */
   class UniformAxialField: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit UniformAxialField();

         /**
          * @brief Simple empty destructor
          */
         virtual ~UniformAxialField();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDFloat amplitude);

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
          * @brief Amplitude of the uniform field
          */
         MHDFloat mA;
   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_UNIFORMAXIALFIELD_HPP
