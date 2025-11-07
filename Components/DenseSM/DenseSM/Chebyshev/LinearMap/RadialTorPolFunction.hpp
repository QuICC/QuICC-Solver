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
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "DenseSM/IGenericProfile.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the generic radial toroidal/poloidal function
 */
class RadialTorPolFunction : public IGenericProfile
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
    * @brief Required spectral truncation
    */
   virtual int nN() const = 0;

   /**
    * @brief Nonzero harmonic degrees
    */
   virtual std::vector<int> ls() const = 0;

   /**
    * @brief Evaluate function on grid
    */
   virtual Internal::Array evaluate(const Internal::Array& r, const int l,
      const int m) const = 0;

   /**
    * @brief Evaluate function derivative on grid
    *
    * @param p Derivative order
    * @param r Radial grid
    * @param l Harmonic degree
    * @param m Harmonic order
    */
   virtual Array evaluateDiff(const int p, const Internal::Array& r, const int l,
      const int m, const Internal::MHDFloat lb, const Internal::MHDFloat ub) const;
protected:
private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_RADIALTORPOLFUNCTION_HPP
