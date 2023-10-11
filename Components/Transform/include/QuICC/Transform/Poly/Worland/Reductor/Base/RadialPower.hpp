/**
 * @file RadialPower.hpp
 * @brief Implementation of the Worland based power spectrum operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_RADIALPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_RADIALPOWER_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandRadialPower.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class RadialPower;

   /**
    * @brief Implementation of the Worland based power spectrum operator on radial grid
    */
   template <>
   class RadialPower<base_t>: public IWorlandRadialPower
   {
      public:
         /**
          * @brief Constructor
          */
         RadialPower();

         /**
          * @brief Destructor
          */
         virtual ~RadialPower() = default;

      protected:
         /**
          * @brief Apply ith operator
          *
          * @param rOut Output radial power
          * @param i    3D mode index
          * @param in   Input spectrum
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

      private:
         /**
          * @brief Make operator
          *
          * @param op         Storage for operator
          * @param igrid      quadrature grid
          * @param iweights   quadrature weights
          * @param i          3D mode index
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;
   };

} // Reductor
} // Worland
} // Poly
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_RADIALPOWER_HPP
