/**
 * @file Power.hpp
 * @brief Implementation of the Worland based power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWER_HPP

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
   class Power;

   /**
    * @brief Implementation of the Worland based power spectrum operator
    */
   template <>
   class Power<base_t>: public IWorlandPower
   {
      public:
         /**
          * @brief Constructor
          */
         Power();

         /**
          * @brief Destructor
          */
         virtual ~Power() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWER_HPP
