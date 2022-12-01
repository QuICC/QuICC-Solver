/**
 * @file PIALegendreOperatorTypes.hpp
 * @brief Interface for a associated Legendre based operator parallel data types
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PIALEGENDREOPERATORTYPES_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PIALEGENDREOPERATORTYPES_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#ifdef QUICC_USE_KOKKOS
#include "KokkosTypedefs.hpp"
#endif

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

/**
 * @brief Interface for a associated Legendre based operator data types
 */
class PIALegendreOperatorTypes {
 public:
   using OpArrayI = QuICC::ArrayI;
   using OpArray = internal::Array;
   using OpMatrix = QuICC::Matrix;
   using OpMatrixZ = QuICC::MatrixZ;
   using OpMatrixR = Eigen::Ref<MatrixZ>;
   using OpMatrixCR = Eigen::Ref<const MatrixZ>;

   using ScalarType = typename OpMatrix::Scalar;

#ifdef QUICC_USE_KOKKOS
   /* #if defined(QUICC_USE_KOKKOS_CUDA) */
   using DataType = Kokkos::complex<double>;
   /* #else
     using DataType = OpMatrixZ::Scalar;
   #endif */

   using OpMatrixLZ = ViewMatrixTypeLeft<DataType>;
   using OpMatrixL = ViewMatrixTypeLeft<ScalarType>;
   using OpMatrixLZH = ViewMatrixTypeLeftHost<OpMatrixZ::Scalar>;

   using OpVectorI = ViewVectorType<int>;
   using OpMatrixI = ViewMatrixTypeLeft<int>;
   using OpVectorITex = ViewVectorTexture<int>;
#endif
};

} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PIALEGENDREOPERATORTYPES_HPP
