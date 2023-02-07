/**
 * @file ScalarYllPerturbation.hpp
 * @brief Scalar field perturbation with Y_l^l spherical harmonic generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SHELL_SCALARYLLPERTURBATION_HPP
#define QUICC_PHYSICAL_KERNEL_SHELL_SCALARYLLPERTURBATION_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

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
          * @param ri         Inner radius
          * @param ro         Outer radius
          * @param epsilon    Amplitude of perturbation
          * @param l          harmonic degree
          */
         void init(const MHDFloat ri, const MHDFloat ro, const MHDFloat epsilon, const int l);

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
          * @brief Amplitude of perturbation
          */
         MHDFloat mEpsilon;

         /**
          * @brief Harmonic degree
          */
         int mL;
   };

} // Shell
} // Kernel
} // Physical
} // QuICC

#endif // QUICC_PHYSICAL_KERNEL_SHELL_SCALARYLLPERTURBATION_HPP
