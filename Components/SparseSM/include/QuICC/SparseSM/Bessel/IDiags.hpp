/**
 * @file IDiags.hpp
 * @brief Interface to spherical Bessel sparse operator diagonals
 */

#ifndef QUICC_SPARSESM_BESSEL_IDIAGS_HPP
#define QUICC_SPARSESM_BESSEL_IDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/BesselKind.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

/**
 * @brief Interface to Bessel sparse operator diagonals
 */
class IDiags
{
public:
   /// Typedef for scalar
   typedef Internal::MHDFloat Scalar_t;

   /// Typedef for coefficient array
   typedef Internal::ACoeff ACoeff_t;

   /**
    * @brief Constructor
    *
    * @param type Type of Bessel basis
    * @param l    Harmonic degree l
    */
   IDiags(const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   virtual ~IDiags() = default;

protected:
   /**
    * @brief Get l
    */
   Scalar_t l() const;

   /**
    * @brief Type of Bessel implementation
    */
   Bessel::BesselKind type() const;

private:
   /**
    * @brief Type of Bessel basis
    */
   BesselKind mType;

   /**
    * @brief l
    */
   Scalar_t mL;
};

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_IDIAGS_HPP
