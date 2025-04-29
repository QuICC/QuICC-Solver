/**
 * @file LlDivS1Dp.hpp
 * @brief Implementation of the associated Legendre based l(l+1)/Sin D_phi integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LLDIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LLDIVS1DP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class LlDivS1Dp;

   /**
    * @brief Implementation of the associated Legendre based l(l+1)/Sin D_phi integrator
    */
   template <>
   class LlDivS1Dp<base_t>: public DivS1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         LlDivS1Dp() = default;

         /**
          * @brief Destructor
          */
         virtual ~LlDivS1Dp() = default;

      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief l(l+1) scaling factors
          */
         virtual void initSpecial() const override;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LLDIVS1DP_HPP
