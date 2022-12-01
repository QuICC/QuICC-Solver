/** 
 * @file Ll2.hpp
 * @brief Implementation of the associated Legendre based l(l+1)^2 P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL2_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL2_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Implementation of the associated Legendre based l(l+1)^2 P integrator
    */ 
   class Ll2: public P<>
   {
      public:
         /**
          * @brief Constructor
          */
         Ll2();

         /**
          * @brief Destructor
          */
         virtual ~Ll2();
         
      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief l(l+1) scaling factors
          */
         virtual void initSpecial() const override;

         /**
          * @brief Storage for l(l+1)^2 factors
          */
         mutable Array mLl2;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL2_HPP
