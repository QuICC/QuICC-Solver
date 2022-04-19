/**
 * @file TorPolHarmonic.hpp
 * @brief Exact Toroidal/Poloidal harmonic state generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_TORPOLHARMONIC_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_TORPOLHARMONIC_HPP

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
    * @brief Exact Toroidal/Poloidal harmonic state generator kernel
    */
   class TorPolHarmonic: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit TorPolHarmonic();

         /**
          * @brief Simple empty destructor
          */
         virtual ~TorPolHarmonic();

         /**
          * @brief Set component modes
          */
         void setModes(const FieldComponents::Spectral::Id compId, const Spectral::Kernel::Complex3DMapType& modes);

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init();

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
          * @brief Shared map of harmonic modes
          */
         Spectral::Kernel::VectSharedZ3DMapType mSHModes;

   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_TORPOLHARMONIC_HPP
