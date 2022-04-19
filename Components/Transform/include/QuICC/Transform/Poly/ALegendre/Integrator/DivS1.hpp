/** 
 * @file DivS1.hpp
 * @brief Implementation of the associated Legendre based 1/Sin integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/IALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Implementation of the associated Legendre based 1/Sin integrator
    */ 
   class DivS1: public IALegendreIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         DivS1();

         /**
          * @brief Destructor
          */
         virtual ~DivS1();
         
      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1_HPP
