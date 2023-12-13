/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "TestSuite/Timestep/RungeKuttaCB/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace RungeKuttaCB {

TestArgs::TestArgs() :
    useDefault(true),
    dumpData(false),
    timeOnly(false),
    ulp(11),
    dim1D(0),
    dim2D(0),
    dim3D(0),
    scheme("")
{}

TestArgs& args()
{
   static TestArgs a;

   return a;
}

} // namespace RungeKuttaCB
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC
