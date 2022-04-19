/**
 * @file BenchmarkTempC1.hpp
 * @brief Temperature benchmark state C1 generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPC1_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPC1_HPP

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
    * @brief Temperature benchmark state C1 generator kernel
    */
   class BenchmarkTempC1: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit BenchmarkTempC1();

         /**
          * @brief Simple empty destructor
          */
         virtual ~BenchmarkTempC1();

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

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPC1_HPP
