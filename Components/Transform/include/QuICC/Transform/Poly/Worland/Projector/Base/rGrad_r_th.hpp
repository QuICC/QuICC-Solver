/**
 * @file rGrad_r_th.hpp
 * @brief Implementation of the Worland based projector of radial derivatives in Grad_r_th
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_RGRAD_R_TH_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_RGRAD_R_TH_HPP

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
   class rGrad_r_th;

   /**
    * @brief Implementation of the Worland based r, th component of the tensor gradient (radial derviatives only)
    */
   template <>
   class rGrad_r_th<base_t>: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         rGrad_r_th();

         /**
          * @brief Destructor
          */
         ~rGrad_r_th() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_BASE_RGRAD_R_TH_HPP
