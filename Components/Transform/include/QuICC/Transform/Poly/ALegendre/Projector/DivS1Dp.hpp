/** 
 * @file DivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/sin P d_phi projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1DP_HPP

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
    * @brief Implementation of the associated Legendre based 1/sin P d_phi projector
    */ 
   class DivS1Dp: public DivS1
   {
      public:
         /**
          * @brief Constructor
          */
         DivS1Dp();

         /**
          * @brief Destructor
          */
         virtual ~DivS1Dp();
         
      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1DP_HPP
