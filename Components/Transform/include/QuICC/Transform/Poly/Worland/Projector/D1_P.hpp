/**
 * @file D1_P.hpp
 * @brief Implementation of the Worland based D projector but 0 mode is P projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1_P_HPP

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
    * @brief Implementation of the Worland based D projector but 0 mode is P projector
    */
   class D1_P: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D1_P();

         /**
          * @brief Destructor
          */
         virtual ~D1_P();

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1_P_HPP
