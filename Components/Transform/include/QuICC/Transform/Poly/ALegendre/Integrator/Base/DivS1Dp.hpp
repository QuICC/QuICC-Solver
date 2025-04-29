/**
 * @file DivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/Sin D_phi integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVS1DP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class DivS1Dp;

   /**
    * @brief Implementation of the associated Legendre based 1/Sin D_phi integrator
    */
   template <>
   class DivS1Dp<base_t>: public DivS1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         DivS1Dp() = default;

         /**
          * @brief Destructor
          */
         virtual ~DivS1Dp() = default;

      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const final;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVS1DP_HPP
