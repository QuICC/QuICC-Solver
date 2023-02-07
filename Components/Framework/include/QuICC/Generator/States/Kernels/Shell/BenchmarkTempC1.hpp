/**
 * @file BenchmarkTempC1.hpp
 * @brief Temperature benchmark state from Christensen's C1 test case generator kernel
 *        DOI: https://doi.org/10.1016/S0031-9201(01)00275-8
 */

#ifndef QUICC_PHYSICAL_KERNEL_SHELL_BENCHMARKTEMPC1_HPP
#define QUICC_PHYSICAL_KERNEL_SHELL_BENCHMARKTEMPC1_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Generator/States/Kernels/Shell/ScalarYllPerturbation.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

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
          * @param ri Inner radius
          * @param ro Outer radius
          */
         void init(const MHDFloat ri, const MHDFloat ro);

      protected:

      private:
   };

} // Shell
} // Kernel
} // Physical
} // QuICC

#endif // QUICC_PHYSICAL_KERNEL_SHELL_BENCHMARKTEMPC1_HPP
