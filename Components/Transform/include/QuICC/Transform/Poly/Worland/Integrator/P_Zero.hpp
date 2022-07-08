/**
 * @file P_Zero.hpp
 * @brief Implementation of the Worland based integrator, zero for l = 0
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_ZERO_HPP

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
    * @brief Implementation of the Worland based integrator, zero for l = 0
    */
   class P_Zero: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         P_Zero();

         /**
          * @brief Destructor
          */
         virtual ~P_Zero();

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_ZERO_HPP
