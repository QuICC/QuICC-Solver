/**
 * @file PIOperatorTypes.hpp
 * @brief Interface for a associated Poly based operator parallel data types
 */

#ifndef QUICC_TRANSFORM_POLY_CUDAIOPERATORTYPES_HPP
#define QUICC_TRANSFORM_POLY_CUDAIOPERATORTYPES_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#ifdef QUICC_USE_KOKKOS
#include "KokkosTypedefs.hpp"
#include "QuICC/Transform/Poly/CudaHipRuntime.hpp"
#endif

namespace QuICC {

namespace Transform {

namespace Poly {

/**
 * @brief Interface for a associated Poly based operator data types
 */
class CudaIOperatorTypes {
 public:

#ifdef QUICC_USE_KOKKOS

   using ScalarType = typename KokkosIOperatorTypes::ScalarType;
   using DataType =  typename KokkosIOperatorTypes::DataType;



   //CUDA TYPES
   //
   // Matrices are stored in column-major order:
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

} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_CUDAIOPERATORTYPES_HPP
