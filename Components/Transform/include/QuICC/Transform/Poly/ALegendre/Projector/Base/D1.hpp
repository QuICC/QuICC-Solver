/**
 * @file D1.hpp
 * @brief Implementation of the associated Legendre based D projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_D1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/IALegendreProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class D1;

   /**
    * @brief Implementation of the associated Legendre based D projector
    */
   template <>
   class D1<base_t>: public IALegendreProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D1() = default;

         /**
          * @brief Destructor
          */
         virtual ~D1() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_D1_HPP
