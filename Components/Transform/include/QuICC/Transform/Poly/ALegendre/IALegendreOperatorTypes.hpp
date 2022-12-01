/**
 * @file IALegendreOperatorTypes.hpp
 * @brief Interface for a associated Legendre based operator data types
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_IALEGENDREOPERATORTYPES_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_IALEGENDREOPERATORTYPES_HPP

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

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

   /**
    * @brief Interface for a associated Legendre based operator data types
    */
   class IALegendreOperatorTypes
   {
      public:
         using OpArrayI = QuICC::ArrayI;
         using OpArray = internal::Array;
         using OpMatrix = QuICC::Matrix;
         using OpMatrixZ = QuICC::MatrixZ;
         using OpMatrixR = Eigen::Ref<MatrixZ>;
         using OpMatrixCR = Eigen::Ref<const MatrixZ>;

         using DataType = typename OpMatrixZ::Scalar;
         using ScalarType = typename OpMatrix::Scalar;

   };

}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_IALEGENDREOPERATORTYPES_HPP
