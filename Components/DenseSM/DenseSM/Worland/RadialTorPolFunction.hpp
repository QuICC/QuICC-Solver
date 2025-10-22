/**
 * @file RadialTorPolFunction.hpp
 * @brief Implementation of the generic radial toroidal/poloidal function
 */

#ifndef QUICC_DENSESM_WORLAND_RADIALTORPOLFUNCTION_HPP
#define QUICC_DENSESM_WORLAND_RADIALTORPOLFUNCTION_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

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

protected:
private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_RADIALTORPOLFUNCTION_HPP
