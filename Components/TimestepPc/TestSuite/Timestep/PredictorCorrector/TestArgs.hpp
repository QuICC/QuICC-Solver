/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_TIMESTEP_PREDICTORCORRECTOR_TESTARGS_HPP
#define QUICC_TESTSUITE_TIMESTEP_PREDICTORCORRECTOR_TESTARGS_HPP

// System includes
//
#include <string>
#include <vector>

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace PredictorCorrector {

struct TestArgs
{
   /// Use default test setup
   bool useDefault;

   /// Write output data to file
   bool dumpData;

   /// Only time execution, don't check data
   bool timeOnly;

   /// Max ulp
   unsigned int ulp;

   /// Dimension 1D
   unsigned int dim1D;

   /// Dimension 1D
   unsigned int dim2D;

   /// Dimension 1D
   unsigned int dim3D;

   /// Timestepping scheme
   std::string scheme;

   /// ID of the tests
   std::vector<int> params;

   /**
    * @brief Constructor
    */
   TestArgs();

   /**
    * @brief Destructor
    */
   ~TestArgs() = default;
};

TestArgs& args();

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_FRAMEWORK_TIMESTEP_PREDICTORCORRECTOR_TESTARGS_HPP
