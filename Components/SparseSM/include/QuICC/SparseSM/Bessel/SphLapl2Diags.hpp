/**
 * @file SphLapl2Diags.hpp
 * @brief Interface to SphLapl2 diagonals for full sphere Bessel SphLapl2 sparse
 * operator
 */

#ifndef QUICC_SPARSESM_BESSEL_SPHLAPL2DIAGS_HPP
#define QUICC_SPARSESM_BESSEL_SPHLAPL2DIAGS_HPP

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
 * @brief Implementation of the full sphere Bessel SphLapl2 sparse operator
 */
class SphLapl2Diags : public IDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param type Type of Bessel basis
    * @param l    Harmonic degree l
    */
   SphLapl2Diags(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   virtual ~SphLapl2Diags() = default;

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

#endif // QUICC_SPARSESM_BESSEL_SPHLAPL2DIAGS_HPP
