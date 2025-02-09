/**
 * @file PowerD1R1.hpp
 * @brief Implementation of the Worland based D R power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWERD1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWERD1R1_HPP

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
   class PowerD1R1;

   /**
    * @brief Implementation of the Worland based D R power spectrum operator
    */
   template <>
   class PowerD1R1<base_t>: public IWorlandPower
   {
      public:
         /**
          * @brief Constructor
          */
         PowerD1R1();

         /**
          * @brief Destructor
          */
         virtual ~PowerD1R1() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_POWERD1R1_HPP
