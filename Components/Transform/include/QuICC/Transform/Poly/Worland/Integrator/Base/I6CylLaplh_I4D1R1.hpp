/**
 * @file I6CylLaplh_I4D1R1.hpp
 * @brief Implementation of the Worland based I6 cylindrical horizontal laplacian integrator but 0 mode is I4 D R1 integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I6CYLLAPLH_I4D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I6CYLLAPLH_I4D1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/IWorlandIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   template <class Impl>
   class I6CylLaplh_I4D1R1;

   /**
    * @brief Implementation of the Worland based I6 cylindrical horizontal laplacian integrator but 0 mode is I4 D R1 integrator
    */
   template <>
   class I6CylLaplh_I4D1R1<base_t>: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I6CylLaplh_I4D1R1() = default;

         /**
          * @brief Destructor
          */
         ~I6CylLaplh_I4D1R1() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_BASE_I6CYLLAPLH_I4D1R1_HPP
