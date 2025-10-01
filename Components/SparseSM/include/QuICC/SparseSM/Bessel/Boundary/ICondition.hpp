/**
 * @file ICondition.hpp
 * @brief Interface to generic Bessel boundary condition
 */

#ifndef QUICC_SPARSESM_BESSEL_BOUNDARY_ICONDITION_HPP
#define QUICC_SPARSESM_BESSEL_BOUNDARY_ICONDITION_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/IDiags.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

/**
 * @brief Interface to generic Bessel boundary condition
 */
class ICondition : public IDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param l       Harmonic degree l
    */
   ICondition(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   virtual ~ICondition() = default;

protected:
private:
};

} // namespace Boundary
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_BOUNDARY_ICONDITION_HPP
