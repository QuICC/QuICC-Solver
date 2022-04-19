/**
 * @file BenchmarkVelC2.hpp
 * @brief Velocity benchmark state C2 generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKVELC2_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKVELC2_HPP

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
    * @brief Velocity benchmark state C2 generator kernel
    */
   class BenchmarkVelC2: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit BenchmarkVelC2();

         /**
          * @brief Simple empty destructor
          */
         virtual ~BenchmarkVelC2();

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

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKVELC2_HPP
