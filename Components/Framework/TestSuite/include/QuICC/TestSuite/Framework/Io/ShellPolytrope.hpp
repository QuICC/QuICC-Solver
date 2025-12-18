/**
 * @file ShellPolytrope.hpp
 * @brief Implementation of the polytropic density field
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_IO_SHELLPOLYTROPE_HPP
#define QUICC_TESTSUITE_FRAMEWORK_IO_SHELLPOLYTROPE_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

/**
 * @brief Implementation of the polytropic field
 */
class ShellPolytrope
    : public QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction
{
public:
   /**
    * @brief Constructor
    */
   ShellPolytrope(const MHDFloat rratio, const MHDFloat Nrho, const MHDFloat npoly);

   /**
    * @brief Destructor
    */
   virtual ~ShellPolytrope() = default;

   /**
    * @brief Spectral truncation
    */
   int nN() const final;

   /**
    * @brief Nonzero harmonic degrees
    */
   std::vector<int> ls() const final;

   /**
    * @brief Evaluate on radial grid
    */
   Internal::Array evaluate(const Internal::Array& r, const int l,
      const int m) const final;

   /**
    * @brief Radius Ratio
    */
   const MHDFloat mRratio;
   /**
    * @brief Scale heights
    */
   const MHDFloat mNrho;
   /**
    * @brief polytropic index
    */
   const MHDFloat mNpoly;


protected:
private:
};

} // namespace Io
} // namespace Framework
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_FRAMEWORK_IO_SHELLPOLYTROPE_HPP
