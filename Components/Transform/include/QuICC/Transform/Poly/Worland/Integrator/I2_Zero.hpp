/** 
 * @file I2_Zero.hpp
 * @brief Implementation of the Worland based I2 integrator but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I2_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I2_ZERO_HPP

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
#include "QuICC/Transform/Poly/Worland/Integrator/IWorlandIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   /**
    * @brief Implementation of the Worland based I2 integrator but 0 mode is zeroed
    */ 
   class I2_Zero: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I2_Zero();

         /**
          * @brief Destructor
          */
         virtual ~I2_Zero();
         
      protected:

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I2_ZERO_HPP
