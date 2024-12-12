/**
 * @file RadialTorPolFunction.hpp
 * @brief Implementation of the generic radial toroidal/poloidal function
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_RADIALTORPOLFUNCTION_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_RADIALTORPOLFUNCTION_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the generic radial toroidal/poloidal function
 */
class RadialTorPolFunction
{
public:
   /**
    * @brief Constructor
    */
   RadialTorPolFunction() = default;

   /**
    * @brief Destructor
    */
   virtual ~RadialTorPolFunction() = default;

   /**
    * @brief Spectral truncation
    */
   virtual int nN() const = 0;

   /**
    * @brief Spectral expansion
    */
   virtual Internal::Array evaluate(const Internal::Array& r, const int l, const int m) const = 0;

protected:

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_RADIALTORPOLFUNCTION_HPP
