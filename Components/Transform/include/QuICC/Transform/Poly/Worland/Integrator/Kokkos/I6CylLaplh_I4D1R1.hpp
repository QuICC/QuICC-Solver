/**
 * @file I6CylLaplh_I4D1R1.hpp
 * @brief Implementation of the Worland based I6 cylindrical horizontal
 * laplacian integrator but 0 mode is I4 D R1 integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I6CYLLAPLH_I4D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I6CYLLAPLH_I4D1R1_HPP


// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/KokkosIWorlandIntegrator.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

template <class Impl> class I6CylLaplh_I4D1R1;

/**
 * @brief Implementation of the Worland based I6 cylindrical horizontal
 * laplacian integrator but 0 mode is I4 D R1 integrator
 */
template <> class I6CylLaplh_I4D1R1<kokkos_t> : public KokkosIWorlandIntegrator
{
public:
   /**
    * @brief Constructor
    */
   I6CylLaplh_I4D1R1() = default;

   /**
    * @brief Destructor
    */
   ~I6CylLaplh_I4D1R1() = default;

   void applyUnitOperator(const OpMatrixLZ& rOut, const OpMatrixLZL& in,
      const OpVectorI& scan, const int totalOpsCols) const final;

protected:
   /**
    * @brief Make operator
    */
   void makeOperator(Matrix& op, const Internal::Array& igrid,
      const Internal::Array& iweights, const int i) const final;
};

} // namespace Integrator
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_PPI6CYLLAPLH_I4D1R1_HPP
