/** 
 * @file LlDivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/sin l(l+1) P d_phi P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1DP_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Projector/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   /**
    * @brief Implementation of the associated Legendre based 1/sin l(l+1) P d_phi projector
    */ 
   class LlDivS1Dp: public DivS1
   {
      public:
         /**
          * @brief Constructor
          */
         LlDivS1Dp();

         /**
          * @brief Destructor
          */
         virtual ~LlDivS1Dp();
         
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
          * @brief Storage for l(l+1) factors
          */
         mutable Array mLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1DP_HPP
