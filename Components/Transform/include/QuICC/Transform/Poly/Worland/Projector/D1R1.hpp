/**
 * @file D1R1.hpp
 * @brief Implementation of the Worland based D R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1R1_HPP

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
    * @brief Implementation of the Worland based D R projector
    */
   class D1R1: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D1R1();

         /**
          * @brief Destructor
          */
         virtual ~D1R1();

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1R1_HPP
