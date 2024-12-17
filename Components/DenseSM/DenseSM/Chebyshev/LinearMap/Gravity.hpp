/**
 * @file Gravity.hpp
 * @brief Implementation of the gravity fields
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_GRAVITY_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_GRAVITY_HPP

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
 * @brief Implementation of the dipolar gravity field
 */
class Gravity: public RadialTorPolFunction
{
public:
   /**
    * @brief Constructor
    */
   Gravity() = default;

   /**
    * @brief Destructor
    */
   virtual ~Gravity() = default;

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

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_GRAVITY_HPP
