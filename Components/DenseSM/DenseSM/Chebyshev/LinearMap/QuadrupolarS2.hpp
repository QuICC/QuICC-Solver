/**
 * @file QuadrupolarS2.hpp
 * @brief Implementation of the quadrupolar S2 field
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_QUADRUPOLARS2_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_QUADRUPOLARS2_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the quadrupolar S2 field
 */
class QuadrupolarS2: public RadialTorPolFunction
{
public:
   /**
    * @brief Constructor
    */
   QuadrupolarS2() = default;

   /**
    * @brief Destructor
    */
   virtual ~QuadrupolarS2() = default;

   /**
    * @brief Spectral truncation
    */
   int nN() const final;

   /**
    * @brief Evaluate on radial grid
    */
   Internal::Array evaluate(const Internal::Array& r, const int l, const int m) const final;

protected:

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_QUADRUPOLARS2_HPP
