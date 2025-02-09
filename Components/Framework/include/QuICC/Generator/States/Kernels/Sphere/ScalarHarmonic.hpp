/**
 * @file ScalarHarmonic.hpp
 * @brief Exact scalar harmonic state generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_SCALARHARMONIC_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_SCALARHARMONIC_HPP

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
#include "QuICC/SpectralKernels/Typedefs.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   /**
    * @brief Exact scalar harmonic state generator kernel
    */
   class ScalarHarmonic: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit ScalarHarmonic();

         /**
          * @brief Simple empty destructor
          */
         ~ScalarHarmonic();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const Spectral::Kernel::Complex3DMapType& shModes);

         /**
          * @brief Compute the physical kernel
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         void compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const final;

      protected:

      private:
         /**
          * @brief Shared list of harmonic modes
          */
         Spectral::Kernel::SharedComplex3DMapType mSHModes;

   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_SCALARHARMONIC_HPP
