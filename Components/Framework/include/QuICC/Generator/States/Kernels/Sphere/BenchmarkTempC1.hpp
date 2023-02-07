/**
 * @file BenchmarkTempC1.hpp
 * @brief Temperature benchmark state C1 generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPC1_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPC1_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Generator/States/Kernels/Sphere/ScalarYllPerturbation.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   /**
    * @brief Temperature benchmark state C1 generator kernel
    */
   class BenchmarkTempC1: public ScalarYllPerturbation
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit BenchmarkTempC1();

         /**
          * @brief Simple empty destructor
          */
         virtual ~BenchmarkTempC1() = default;

         /**
          * @brief Initialize kernel
          *
          * @param amplitude_bg Amplitude of background state
          * @param epsilon      Amplitude of perturbation
          */
         void init(const MHDFloat amplitude_bg, const MHDFloat epsilon);

      protected:

      private:
   };

} // Sphere
} // Kernel
} // Physical
} // QuICC

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_BENCHMARKTEMPC1_HPP
