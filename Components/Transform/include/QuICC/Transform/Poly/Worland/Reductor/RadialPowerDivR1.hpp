/** 
 * @file RadialPowerDivR1.hpp
 * @brief Implementation of the Worland based 1/R power spectrum operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandRadialPower.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based 1/R power spectrum operator on radial grid
    */ 
   class RadialPowerDivR1: public IWorlandRadialPower
   {
      public:
         /**
          * @brief Constructor
          */
         RadialPowerDivR1();

         /**
          * @brief Destructor
          */
         virtual ~RadialPowerDivR1() = default;
         
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
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const override;
   };

} // Reductor
} // Worland
} // Poly
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1_HPP
