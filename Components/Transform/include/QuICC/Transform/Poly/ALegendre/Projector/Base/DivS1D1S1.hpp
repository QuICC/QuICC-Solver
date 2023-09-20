/**
 * @file DivS1D1S1.hpp
 * @brief Implementation of the associated Legendre based 1/sin D sin projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_DivS1D1S1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_DivS1D1S1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/IALegendreProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class DivS1D1S1;

   /**
    * @brief Implementation of the associated Legendre based 1/sin D sin projector
    */
   template <>
   class DivS1D1S1<base_t>: public IALegendreProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivS1D1S1() = default;

         /**
          * @brief Destructor
          */
         virtual ~DivS1D1S1() = default;

      protected:

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_DivS1D1S1_HPP
