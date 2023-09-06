/**
 * @file ScalarYllPerturbation.hpp
 * @brief Scalar field perturbation with Y_l^l spherical harmonic generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_SCALARYLLPERTURBATION_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_SCALARYLLPERTURBATION_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   /**
    * @brief Scalar field perturbation with m-fold symmetry generator kernel
    */
   class ScalarYllPerturbation: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit ScalarYllPerturbation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ScalarYllPerturbation() = default;

         /**
          * @brief Initialize kernel
          *
          * @param amplitude_bg Amplitude of background state
          * @param epsilon      Amplitude of perturbation
          * @param l          harmonic degree
          */
         void init(const MHDFloat amplitude_bg, const MHDFloat epsilon, const int l);

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
          * @brief Amplitude of background state
          */
         MHDFloat mBg;

         /**
          * @brief Amplitude of perturbation
          */
         MHDFloat mEpsilon;

         /**
          * @brief Harmonic degree
          */
         int mL;
   };

} // Sphere
} // Kernel
} // Physical
} // QuICC

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_SCALARYLLPERTURBATION_HPP
