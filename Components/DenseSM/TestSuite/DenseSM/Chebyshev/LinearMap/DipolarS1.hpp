/**
 * @file DipolarS1.hpp
 * @brief Implementation of the dipolar S1 field
 */

#ifndef QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_DIPOLARS1_HPP
#define QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_DIPOLARS1_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the dipolar S1 field
 */
class DipolarS1: public QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction
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

}
}
}
}
} // namespace QuICC

#endif // QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_DIPOLARS1_HPP
