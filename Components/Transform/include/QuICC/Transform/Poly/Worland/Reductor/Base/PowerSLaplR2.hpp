/**
 * @file PowerSLaplR2.hpp
 * @brief Implementation of the Worland based spherical Laplacian R^2 power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWERSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWERSLAPLR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandPower.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class PowerSLaplR2;

   /**
    * @brief Implementation of the Worland based Spherical Laplacian R^2 power spectrum operator
    */
   template <>
   class PowerSLaplR2<base_t>: public IWorlandPower
   {
      public:
         /**
          * @brief Constructor
          */
         PowerSLaplR2();

         /**
          * @brief Destructor
          */
         virtual ~PowerSLaplR2() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const int i) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWERSLAPLR2_HPP
