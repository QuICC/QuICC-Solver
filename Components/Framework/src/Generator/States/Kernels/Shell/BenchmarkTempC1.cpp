/**
 * @file BenchmarkTempC1.cpp
 * @brief Source of benchmark state C1 for temperature generator kernel
 */

// System includes
//

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Generator/States/Kernels/Shell/BenchmarkTempC1.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

   BenchmarkTempC1::BenchmarkTempC1()
      : ScalarYllPerturbation()
   {
   }

   void BenchmarkTempC1::init(const MHDFloat ri, const MHDFloat ro)
   {
      ScalarYllPerturbation::init(ri, ro, 1./5., 4);
   }

} // Shell
} // Kernel
} // Physical
} // QuICC
