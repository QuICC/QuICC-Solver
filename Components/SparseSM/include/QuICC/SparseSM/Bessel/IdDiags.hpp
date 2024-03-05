/**
 * @file IdDiags.hpp
 * @brief Interface to I2 diagonals for full sphere Bessel (restricted) identity
 * sparse operator
 */

#ifndef QUICC_SPARSESM_BESSEL_IDDIAGS_HPP
#define QUICC_SPARSESM_BESSEL_IDDIAGS_HPP

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
 * @brief Implementation of the full sphere Bessel (restricted) identity sparse
 * operator
 */
class IdDiags : public IDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param type Type of Bessel basis
    * @param l    Harmonic degree l
    */
   IdDiags(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   virtual ~IdDiags() = default;

   /**
    * @brief Main diagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d0(const ACoeff_t& n) const;

protected:
private:
};

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_IDDIAGS_HPP
