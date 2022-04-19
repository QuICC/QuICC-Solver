/** 
 * @file DivLlD1.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) D integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLD1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLD1_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/D1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Implementation of the associated Legendre based 1/(l(l+1)) D integrator
    */ 
   class DivLlD1: public D1
   {
      public:
         /**
          * @brief Constructor
          */
         DivLlD1();

         /**
          * @brief Destructor
          */
         virtual ~DivLlD1();
         
      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;

         /**
          * @brief l(l+1) scaling factors
          */
         virtual void initSpecial() const;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLD1_HPP
