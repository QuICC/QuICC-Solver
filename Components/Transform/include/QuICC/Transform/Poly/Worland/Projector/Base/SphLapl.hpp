/**
 * @file SphLapl.hpp
 * @brief Implementation of the Worland based spherical laplacian projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_SPHLAPL_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_SPHLAPL_HPP

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
   class SphLapl;

   /**
    * @brief Implementation of the Worland based spherical laplacian projector
    */
   template <>
   class SphLapl<base_t>: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         SphLapl();

         /**
          * @brief Destructor
          */
         ~SphLapl() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_SPHLAPL_HPP
