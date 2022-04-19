/**
 * @file SolidBodyRotation.hpp
 * @brief Solid body rotation generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_SOLIDBODYROTATION_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_SOLIDBODYROTATION_HPP

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
    * @brief Solid body rotating generator kernel
    */
   class SolidBodyRotation: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit SolidBodyRotation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~SolidBodyRotation();

         /**
          * @brief Initialize kernel
          */
         void init(const MHDFloat x, const MHDFloat y, const MHDFloat z);

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
          * @brief X component
          */
         MHDFloat mX;

         /**
          * @brief Y component
          */
         MHDFloat mY;

         /**
          * @brief Z component
          */
         MHDFloat mZ;
   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_SOLIDBODYROTATION_HPP
