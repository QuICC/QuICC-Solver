/**
 * @file DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based 1/R D R projector and zero for l = 0
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_ZERO_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/IWorlandProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   /**
    * @brief Implementation of the Worland based 1/R D R projector and zero for l = 0
    */
   class DivR1D1R1_Zero: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1D1R1_Zero();

         /**
          * @brief Destructor
          */
         virtual ~DivR1D1R1_Zero();

      protected:

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_ZERO_HPP
