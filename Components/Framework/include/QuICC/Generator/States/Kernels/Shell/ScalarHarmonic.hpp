/**
 * @file ScalarHarmonic.hpp
 * @brief Exact scalar harmonic state generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SHELL_SCALARHARMONIC_HPP
#define QUICC_PHYSICAL_KERNEL_SHELL_SCALARHARMONIC_HPP

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

namespace Shell {

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
         virtual ~ScalarHarmonic();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDFloat ri, const MHDFloat ro, const Spectral::Kernel::Complex3DMapType& shModes);

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
          * @brief Inner radius of spherical shell
          */
         MHDFloat mRi;

         /**
          * @brief Outer radius of spherical shell
          */
         MHDFloat mRo;

         /**
          * @brief Shared list of harmonic modes
          */
         Spectral::Kernel::SharedComplex3DMapType mSHModes;

   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SHELL_SCALARHARMONIC_HPP
