/**
 * @file CylLaplh.hpp
 * @brief Implementation of the Worland based cylindrical horizontal laplacian projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_CYLLAPLH_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_CYLLAPLH_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/IWorlandProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   template <class Impl>
   class CylLaplh;

   /**
    * @brief Implementation of the Worland based cylindrical horizontal laplacian projector
    */
   template <>
   class CylLaplh<base_t>: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         CylLaplh() = default;

         /**
          * @brief Destructor
          */
         ~CylLaplh() = default;

      protected:

      private:
         /**
          * @brief Make operator
          */
         void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const final;

         /**
          * @brief Apply ith operator
          */
         void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const final;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_CYLLAPLH_HPP
