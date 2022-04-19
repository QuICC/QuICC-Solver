/** 
 * @file I4DivR1D1R1_I2.hpp
 * @brief Implementation of the Worland based I4 1/R1 D R1 integrator but 0 mode is I2 P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1D1R1_I2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1D1R1_I2_HPP

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
    * @brief Implementation of the Worland based I4 1/R1 D R1 integrator but 0 mode is I2 P integrator
    */ 
   class I4DivR1D1R1_I2: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I4DivR1D1R1_I2();

         /**
          * @brief Destructor
          */
         virtual ~I4DivR1D1R1_I2();
         
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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1D1R1_I2_HPP
