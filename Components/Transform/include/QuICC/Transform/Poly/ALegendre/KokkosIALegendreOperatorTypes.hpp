/**
 * @file KokkosIALegendreOperatorTypes.hpp
 * @brief Interface for a associated Legendre based operator parallel data types
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_KOKKOSIALEGENDREOPERATORTYPES_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_KOKKOSIALEGENDREOPERATORTYPES_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "KokkosTypedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

/**
 * @brief Interface for a associated Legendre based operator data types
 */
class KokkosIALegendreOperatorTypes {
 public:
   using ScalarType = QuICC::Matrix::Scalar;

   /* #if defined(QUICC_USE_KOKKOS_CUDA) */
   using DataType = Kokkos::complex<double>;
   /* #else
     using DataType = QuICC::MatrixZ::Scalar;
   #endif */

   using OpMatrixLZ = ViewMatrixTypeLeft<DataType>;
   using OpMatrixL = ViewMatrixTypeLeft<ScalarType>;
   using OpMatrixLZH = ViewMatrixTypeLeftHost<QuICC::MatrixZ::Scalar>;

   using OpVectorI = ViewVectorType<int>;
   using OpMatrixI = ViewMatrixTypeLeft<int>;

};

} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_KOKKOSIALEGENDREOPERATORTYPES_HPP
