/**
 * @file TestUtils.hpp
 * @brief Test model setup to validate IO
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_IO_TESTUTILS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_IO_TESTUTILS_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

struct TestParameters
{
   std::shared_ptr<Resolution> spRes;
   MHDFloat maxUlp;
   std::string datadir;
   std::string refdir;
   std::vector<::QuICC::Io::Variable::SharedIVariableAsciiWriter> files;
};

template <typename TScheme> void testFile(TestParameters& test);

// Typedef for storing result of error checks
typedef std::tuple<bool,MHDFloat,MHDFloat> ErrorType;

/// Create the spatial scheme
template <typename TScheme> std::shared_ptr<TScheme> createScheme();

/// create the runner
std::shared_ptr<StateGenerator> createRunner(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme, TestParameters& test);

/// Create reference states
QuICC::Equations::SharedIEquation createStates(std::shared_ptr<StateGenerator> spRunner);

/// Check files
void checkFiles(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme, const TestParameters& test);

/// Check sphere files
std::vector<std::tuple<std::string,int,int,int>> checkSphereFiles(const TestParameters& test);

/// Check sphere files
std::vector<std::tuple<std::string,int,int,int>> checkShellFiles(const TestParameters& test);

/**
 * @brief Compute ULP
 */
ErrorType computeUlp(const MHDFloat data, const MHDFloat ref, const MHDFloat refMod, const MHDFloat tol, const MHDFloat eps);

template <typename TScheme> std::shared_ptr<TScheme> createScheme()
{
   // Set and tune spatial scheme
   auto spScheme = std::make_shared<TScheme>(VectorFormulation::TORPOL, GridPurpose::SIMULATION);

   return spScheme;
}

template <typename TScheme> void testFile(TestParameters& test)
{
   // Register IDs
   QuICC::PhysicalNames::registerAll();

   int status = 0;

   auto spScheme = createScheme<TScheme>();

   // Create simulation
   decltype(createRunner(spScheme, test)) spRunner;

   // Exception handling during the initialisation part
   try
   {
      // Create state generator
      spRunner = createRunner(spScheme, test);
   }

   // If exception is thrown, finalise (close files) and return
   catch(std::logic_error& e)
   {
      try
      {
         QuICC::QuICCEnv().abort(e.what());
      }
      catch(std::logic_error& ee)
      {
         std::cerr << ee.what() << std::endl;
      }

      status = -1;
   }

   if(status == 0)
   {
      // Run the simulation
      spRunner->run();
   }

   // Check output
   checkFiles(spScheme, test);

   // Cleanup and close file handles
   spRunner->finalize();
}

} // namespace Io
} // namespace Framework
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_FRAMEWORK_IO_TESTUTILS_HPP
