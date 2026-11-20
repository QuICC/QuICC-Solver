/**
 * @file BenchmarkVelSteady.hpp
 * @brief Velocity benchmark state for steady dynamo solution generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKVELSTEADY_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKVELSTEADY_HPP

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
    * @brief Velocity benchmark state for steady dynamo solution  generator kernel
    */
   class BenchmarkVelSteady: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit BenchmarkVelSteady() = default;

         /**
          * @brief Simple empty destructor
          */
         virtual ~BenchmarkVelSteady() = default;

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
   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKVELSTEADY_HPP
