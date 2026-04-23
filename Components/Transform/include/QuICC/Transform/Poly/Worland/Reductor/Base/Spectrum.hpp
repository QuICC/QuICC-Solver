/**
 * @file Spectrum.hpp
 * @brief Implementation of the Worland based spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_Spectrum_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_Spectrum_HPP

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
   class Spectrum;

   /**
    * @brief Implementation of the Worland based R^2 power spectrum operator
    */
   template <>
   class Spectrum<base_t>: public IWorlandPower
   {
      public:
         /**
          * @brief Constructor
          */
         Spectrum();

         /**
          * @brief Destructor
          */
         virtual ~Spectrum() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_Spectrum_HPP
