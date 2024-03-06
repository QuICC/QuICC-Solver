/**
 * @file KokkosIOperatorTypes.hpp
 * @brief Interface for a quadrature based operator parallel data types
 */

#ifndef QUICC_TRANSFORM_POLY_KOKKOSIOPERATORTYPES_HPP
#define QUICC_TRANSFORM_POLY_KOKKOSIOPERATORTYPES_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "KokkosTypedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

/**
 * @brief Interface for a associated based operator data types
 */
class KokkosIOperatorTypes {
 public:
   using ScalarType = Matrix::Scalar;

   using DataType = Kokkos::complex<double>;

   using OpMatrixLZ = ViewMatrixType<DataType>;
   using OpMatrixLZL = ViewMatrixType<DataType>;
   using OpMatrixL = ViewMatrixType<ScalarType>;
   using OpMatrixLZH = ViewMatrixTypeHost<MatrixZ::Scalar>;

   using OpVectorI = ViewVectorType<int>;
   using OpMatrixI = ViewMatrixType<int>;

};

} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_KOKKOSIOPERATORTYPES_HPP
