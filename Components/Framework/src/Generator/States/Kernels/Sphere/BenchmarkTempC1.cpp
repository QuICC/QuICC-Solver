/**
 * @file BenchmarkTempC1.cpp
 * @brief Source of benchmark state C1 for temperature generator kernel
 */

// System includes
//

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempC1.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   BenchmarkTempC1::BenchmarkTempC1()
      : ScalarYllPerturbation()
   {
   }

   void BenchmarkTempC1::init(const MHDFloat amplitude_bg, const MHDFloat epsilon)
   {
      ScalarYllPerturbation::init(amplitude_bg, epsilon, 3);
   }

} // Sphere
} // Kernel
} // Physical
} // QuICC
