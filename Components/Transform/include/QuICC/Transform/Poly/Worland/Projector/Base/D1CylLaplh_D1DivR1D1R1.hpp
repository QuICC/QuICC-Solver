/**
 * @file D1CylLaplh_D1DivR1D1R1.hpp
 * @brief Implementation of the Worland based D of cylindrical horizontal laplacian projector but 0 mode is D 1/R D R
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_D1CYLLAPLH_D1DIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_D1CYLLAPLH_D1DIVR1D1R1_HPP

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
   class D1CylLaplh_D1DivR1D1R1;

   /**
    * @brief Implementation of the Worland based D of cylindrical horizontal laplacian projector but 0 mode is D 1/R D R
    */
   template <>
   class D1CylLaplh_D1DivR1D1R1<base_t>: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D1CylLaplh_D1DivR1D1R1() = default;

         /**
          * @brief Destructor
          */
         ~D1CylLaplh_D1DivR1D1R1() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_D1CYLLAPLH_D1DIVR1D1R1_HPP
