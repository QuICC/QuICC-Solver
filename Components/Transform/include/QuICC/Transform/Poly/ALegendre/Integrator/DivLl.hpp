/** 
 * @file DivLl.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLL_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLL_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Implementation of the associated Legendre based 1/(l(l+1)) P integrator
    */ 
   class DivLl: public P<>
   {
      public:
         /**
          * @brief Constructor
          */
         DivLl();

         /**
          * @brief Destructor
          */
         virtual ~DivLl();
         
      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const final;

         /**
          * @brief l(l+1) scaling factors
          */
         void initSpecial() const;

         /**
          * @brief Storage for 1/l(l+1) factors
          */
         mutable Array mDivLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLL_HPP
