/** 
 * @file Power.hpp
 * @brief Implementation of the Worland based power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWER_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandPower.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based power spectrum operator
    */ 
   class Power: public IWorlandPower
   {
      public:
         /**
          * @brief Constructor
          */
         Power();

         /**
          * @brief Destructor
          */
         virtual ~Power();
         
      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const int i) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWER_HPP
