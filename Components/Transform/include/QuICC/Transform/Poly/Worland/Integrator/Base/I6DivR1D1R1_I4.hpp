/**
 * @file I6DivR1D1R1_I4.hpp
 * @brief Implementation of the Worland based I6 1/R1 D R1 integrator but 0 mode is I4 P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I6DIVR1D1R1_I4_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I6DIVR1D1R1_I4_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/IWorlandIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   template <class Impl>
   class I6DivR1D1R1_I4;

   /**
    * @brief Implementation of the Worland based I6 1/R1 D R1 integrator but 0 mode is I4 P integrator
    */
   template <>
   class I6DivR1D1R1_I4<base_t>: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I6DivR1D1R1_I4() = default;

         /**
          * @brief Destructor
          */
         ~I6DivR1D1R1_I4() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE _I6DIVR1D1R1_I4_HPP
