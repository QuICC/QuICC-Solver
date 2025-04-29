/**
 * @file SphLaplDiags.hpp
 * @brief Interface to SphLapl diagonals for full sphere Bessel SphLapl sparse
 * operator
 */

#ifndef QUICC_SPARSESM_BESSEL_INSULATING_SPHLAPLDIAGS_HPP
#define QUICC_SPARSESM_BESSEL_INSULATING_SPHLAPLDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/SphLaplDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Insulating {

/**
 * @brief Implementation of the full sphere Bessel SphLapl sparse operator
 */
class SphLaplDiags : public QuICC::SparseSM::Bessel::SphLaplDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param l Harmonic degree
    */
   SphLaplDiags(const int l);

   /**
    * @brief Destructor
    */
   virtual ~SphLaplDiags() = default;

   /**
    * @brief Main diagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d0(const ACoeff_t& n) const final;

protected:
private:
};

} // namespace Insulating
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_INSULATING_SPHLAPLDIAGS_HPP
