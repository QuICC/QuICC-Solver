/**
 * @file PIALegendreOperatorTypes.hpp
 * @brief Interface for a associated Legendre based operator parallel data types
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAIALEGENDREOPERATORTYPES_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAIALEGENDREOPERATORTYPES_HPP

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

#ifdef QUICC_USE_KOKKOS_CUDA
#include <cuComplex.h>
#endif

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

/**
 * @brief Interface for a associated Legendre based operator data types
 */
class CudaIALegendreOperatorTypes {
 public:
   using OpArrayI = QuICC::ArrayI;
   using OpArray = internal::Array;
   using OpMatrix = QuICC::Matrix;
   using OpMatrixZ = QuICC::MatrixZ;
   using OpMatrixR = Eigen::Ref<MatrixZ>;
   using OpMatrixCR = Eigen::Ref<const MatrixZ>;

   using ScalarType = typename OpMatrix::Scalar;

#ifdef QUICC_USE_KOKKOS_CUDA
   /* using DataType = cuDoubleComplex; */
#endif

#ifdef QUICC_USE_KOKKOS

   using DataType = Kokkos::complex<double>;

   using OpMatrixLZ = ViewMatrixTypeLeft<DataType>;
   using OpMatrixL = ViewMatrixTypeLeft<ScalarType>;
   using OpMatrixLZH = ViewMatrixTypeLeftHost<OpMatrixZ::Scalar>;

   using OpVectorI = ViewVectorType<int>;
   using OpMatrixI = ViewMatrixTypeLeft<int>;
   using OpVectorITex = ViewVectorTexture<int>;

   template<typename T = float>
   struct CuMatrix {
      int width;
      int height;
      int stride;
      int block_stride;
      T *elements;
      using Scalar = T;
   };

   using OpMatrixZC = CuMatrix<DataType>;
   using OpMatrixC = CuMatrix<ScalarType>;

#endif
};

} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PIALEGENDREOPERATORTYPES_HPP
