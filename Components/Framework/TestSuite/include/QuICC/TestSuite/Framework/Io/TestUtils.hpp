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

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

struct TestParameters
{
   std::shared_ptr<Resolution> spRes;
   MHDFloat maxUlp;
};

// Typedef for storing result of error checks
typedef std::tuple<bool,MHDFloat,MHDFloat> ErrorType;

/// Create the spatial scheme
template <typename TScheme> std::shared_ptr<TScheme> createScheme();

/// create the runner
std::shared_ptr<StateGenerator> createRunner(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme, TestParameters& test);

/// Create reference states
QuICC::Equations::SharedIEquation createStates(std::shared_ptr<StateGenerator> spRunner);

/// Add ASCII files for sphere schemes
void addSphereFiles(std::shared_ptr<StateGenerator> spRunner);

/// Add ASCII files for shell schemes
void addShellFiles(std::shared_ptr<StateGenerator> spRunner);

/// Check files
void checkFiles(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme, const TestParameters& test);

/// Check sphere files
void checkSphereFiles(const TestParameters& test);

/// Check sphere files
void checkShellFiles(const TestParameters& test);

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

} // namespace Io
} // namespace Framework
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_FRAMEWORK_IO_TESTUTILS_HPP
