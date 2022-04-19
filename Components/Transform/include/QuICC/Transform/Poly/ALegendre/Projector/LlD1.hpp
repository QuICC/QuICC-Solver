/** 
 * @file LlD1.hpp
 * @brief Implementation of the associated Legendre based l(l+1) D projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLD1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLD1_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Projector/D1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   /**
    * @brief Implementation of the associated Legendre based l(l+1) D projector
    */ 
   class LlD1: public D1
   {
      public:
         /**
          * @brief Constructor
          */
         LlD1();

         /**
          * @brief Destructor
          */
         virtual ~LlD1();
         
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
          * @brief Storage for l dependent factors
          */
         mutable Array mLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLD1_HPP
