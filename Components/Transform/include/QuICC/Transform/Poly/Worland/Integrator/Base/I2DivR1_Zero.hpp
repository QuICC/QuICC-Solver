/**
 * @file I2DivR1_Zero.hpp
 * @brief Implementation of the Worland based I2 1/R1 integrator but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I2DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I2DIVR1_ZERO_HPP

// System includes
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
   class I2DivR1_Zero;

   /**
    * @brief Implementation of the Worland based I2 1/R1 integrator but 0 mode is zeroed
    */
   template <>
   class I2DivR1_Zero<base_t>: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I2DivR1_Zero();

         /**
          * @brief Destructor
          */
         ~I2DivR1_Zero() = default;

      protected:

      private:
         /**
          * @brief Make operator
          */
         void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const final;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I2DIVR1_ZERO_HPP
