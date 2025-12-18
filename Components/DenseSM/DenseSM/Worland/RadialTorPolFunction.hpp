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
#include "DenseSM/IGenericProfile.hpp"


namespace QuICC {

namespace DenseSM {

namespace Worland {

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

protected:
private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_RADIALTORPOLFUNCTION_HPP
