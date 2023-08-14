/**
 * @file P.hpp
 * @brief Implementation of the Worland based integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_P_HPP

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
   class P;

   /**
    * @brief Implementation of the Worland based integrator
    */
   template <>
   class P<base_t>: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         P();

         /**
          * @brief Destructor
          */
         ~P() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_P_HPP
