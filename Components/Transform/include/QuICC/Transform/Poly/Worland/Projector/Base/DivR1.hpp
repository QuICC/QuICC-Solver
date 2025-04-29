/**
 * @file DivR1.hpp
 * @brief Implementation of the Worland based 1/R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_DIVR1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_DIVR1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/IWorlandProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   template <class Impl>
   class DivR1;

   /**
    * @brief Implementation of the Worland based 1/R projector
    */
   template <>
   class DivR1<base_t>: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1();

         /**
          * @brief Destructor
          */
         ~DivR1() = default;

      protected:

      private:
         /**
          * @brief Make operator
          */
         void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const final;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_DIVR1_HPP
