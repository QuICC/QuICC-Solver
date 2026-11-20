/**
 * @file BenchmarkTempSteady.hpp
 * @brief Temperature benchmark state for steady dynamo solution generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPSTEADY_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPSTEADY_HPP

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
    * @brief Temperature benchmark state for steady dynamo solution generator kernel
    */
   class BenchmarkTempSteady: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit BenchmarkTempSteady() = default;

         /**
          * @brief Simple empty destructor
          */
         virtual ~BenchmarkTempSteady() = default;

         /**
          * @brief Initialize kernel
          *
          * @param amplitude_bg Amplitude of background state
          * @param epsilon      Amplitude of perturbation
          */
         void init(const MHDFloat amplitude_bg, const MHDFloat epsilon);

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
   };

} // Sphere
} // Kernel
} // Physical
} // QuICC

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPSTEADY_HPP
