/**
 * @file P_Zero.hpp
 * @brief Implementation of the Worland based integrator, zero for l = 0
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_P_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_P_ZERO_HPP

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
   class P_Zero;

   /**
    * @brief Implementation of the Worland based integrator, zero for l = 0
    */
   template <>
   class P_Zero<base_t>: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         P_Zero();

         /**
          * @brief Destructor
          */
         ~P_Zero() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_P_ZERO_HPP
