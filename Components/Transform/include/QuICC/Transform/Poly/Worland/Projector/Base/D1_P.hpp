/**
 * @file D1_P.hpp
 * @brief Implementation of the Worland based D projector but 0 mode is P projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_D1_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_D1_P_HPP

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
   class D1_P;

   /**
    * @brief Implementation of the Worland based D projector but 0 mode is P projector
    */
   template <>
   class D1_P<base_t>: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D1_P() = default;

         /**
          * @brief Destructor
          */
         ~D1_P() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_D1_P_HPP
