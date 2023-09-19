/**
 * @file CudaIALegendreOperatorTypes.hpp
 * @brief Interface for a associated Legendre based operator parallel data types
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAIALEGENDREOPERATORTYPES_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAIALEGENDREOPERATORTYPES_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "KokkosTypedefs.hpp"

#include <cuComplex.h>

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

/**
 * @brief Interface for a associated Legendre based operator data types
 */
class CudaIALegendreOperatorTypes {
public:

   using ScalarType = typename KokkosIALegendreOperatorTypes::ScalarType;
   using DataType = typename KokkosIALegendreOperatorTypes::DataType;

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

};

} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PIALEGENDREOPERATORTYPES_HPP
