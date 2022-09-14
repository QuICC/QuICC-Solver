/** 
 * @file RadialPower.hpp
 * @brief Implementation of the Worland based power spectrum operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWER_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandRadialPower.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based power spectrum operator on radial grid
    */ 
   class RadialPower: public IWorlandRadialPower
   {
      public:
         /**
          * @brief Constructor
          */
         RadialPower();

         /**
          * @brief Destructor
          */
         virtual ~RadialPower();
         
      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWER_HPP
