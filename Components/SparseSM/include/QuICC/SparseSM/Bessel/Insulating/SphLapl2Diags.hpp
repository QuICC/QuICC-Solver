/**
 * @file SphLapl2Diags.hpp
 * @brief Interface to SphLapl2 diagonals for full sphere Bessel SphLapl2 sparse
 * operator
 */

#ifndef QUICC_SPARSESM_BESSEL_INSULATING_SPHLAPL2DIAGS_HPP
#define QUICC_SPARSESM_BESSEL_INSULATING_SPHLAPL2DIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/SphLapl2Diags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Insulating {

/**
 * @brief Implementation of the full sphere Bessel SphLapl2 sparse operator
 */
class SphLapl2Diags : public QuICC::SparseSM::Bessel::SphLapl2Diags
{
public:
   /**
    * @brief Constructor
    *
    * @param l Harmonic degree
    */
   SphLapl2Diags(const int l);

   /**
    * @brief Destructor
    */
   virtual ~SphLapl2Diags() = default;

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

#endif // QUICC_SPARSESM_BESSEL_INSULATING_SPHLAPL2DIAGS_HPP
