/**
 * @file DivS1.hpp
 * @brief Implementation of the associated Legendre based 1/Sin integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVS1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/IALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class DivS1;

   /**
    * @brief Implementation of the associated Legendre based 1/Sin integrator
    */
   template <>
   class DivS1<base_t>: public IALegendreIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         DivS1() = default;

         /**
          * @brief Destructor
          */
         virtual ~DivS1() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

       private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVS1_HPP
