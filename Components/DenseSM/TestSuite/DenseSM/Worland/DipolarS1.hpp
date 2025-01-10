/**
 * @file DipolarS1.hpp
 * @brief Implementation of the dipolar S1 field
 */

#ifndef QUICC_TESTSUITE_DENSESM_WORLAND_DIPOLARS1_HPP
#define QUICC_TESTSUITE_DENSESM_WORLAND_DIPOLARS1_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the dipolar S1 field
 */
class DipolarS1: public QuICC::DenseSM::Worland::RadialTorPolFunction
{
public:
   /**
    * @brief Constructor
    */
   DipolarS1() = default;

   /**
    * @brief Destructor
    */
   virtual ~DipolarS1() = default;

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
   Internal::Array evaluate(const Internal::Array& r, const int l, const int m) const final;

protected:

private:
};

} // Worland
} // DenseSM
} // TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_DENSESM_WORLAND_DIPOLARS1_HPP
