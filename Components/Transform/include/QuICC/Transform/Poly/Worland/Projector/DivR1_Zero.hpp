/**
 * @file DivR1_Zero.hpp
 * @brief Implementation of the Worland based 1/R projector but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1_ZERO_HPP

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
    * @brief Implementation of the Worland based 1/R projector but 0 mode is zeroed
    */
   class DivR1_Zero: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1_Zero();

         /**
          * @brief Destructor
          */
         virtual ~DivR1_Zero();

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1_ZERO_HPP
