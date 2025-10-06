/**
 * @file Polytrope.hpp
 * @brief Implementation of the polytropic density field
 */

#ifndef QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_POLYTROPE_HPP
#define QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_POLYTROPE_HPP

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
 * @brief Implementation of the polytropic field
 */
class Polytrope: public QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction
{
public:
   /**
    * @brief Constructor
    */
   Polytrope(const MHDFloat rratio, const MHDFloat Nrho, const MHDFloat npoly);

   /**
    * @brief Destructor
    */
   virtual ~Polytrope() = default;

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

}
}
}
}
} // namespace QuICC

#endif // QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_POLYTROPE_HPP
