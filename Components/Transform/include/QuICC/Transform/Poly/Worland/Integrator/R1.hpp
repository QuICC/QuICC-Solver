/** 
 * @file R1.hpp
 * @brief Implementation of the Worland based R1 integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_R1_HPP

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
    * @brief Implementation of the Worland based R1 integrator
    */ 
   class R1: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         R1();

         /**
          * @brief Destructor
          */
         virtual ~R1();
         
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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_R1_HPP
