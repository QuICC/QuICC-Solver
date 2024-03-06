/**
 * @file SphLaplDiags.hpp
 * @brief Interface to SphLapl diagonals for full sphere Bessel SphLapl sparse
 * operator
 */

#ifndef QUICC_SPARSESM_BESSEL_SPHLAPLDIAGS_HPP
#define QUICC_SPARSESM_BESSEL_SPHLAPLDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/IDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

/**
 * @brief Implementation of the full sphere Bessel SphLapl sparse operator
 */
class SphLaplDiags : public IDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param type Type of Bessel basis
    * @param l    Harmonic degree l
    */
   SphLaplDiags(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   virtual ~SphLaplDiags() = default;

   /**
    * @brief Main diagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d0(const ACoeff_t& n) const = 0;

protected:
private:
};

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_SPHLAPLDIAGS_HPP
