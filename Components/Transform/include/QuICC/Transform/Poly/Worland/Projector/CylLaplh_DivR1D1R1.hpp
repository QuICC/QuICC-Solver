/**
 * @file CylLaplh_DivR1D1R1.hpp
 * @brief Implementation of the Worland based cylindrical horizontal laplacian projector but 0 mode 1/R D R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_DIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_DIVR1D1R1_HPP

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
    * @brief Implementation of the Worland based cylindrical horizontal laplacian projector but 0 mode 1/R D R projector
    */
   class CylLaplh_DivR1D1R1: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         CylLaplh_DivR1D1R1();

         /**
          * @brief Destructor
          */
         virtual ~CylLaplh_DivR1D1R1();

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_DIVR1D1R1_HPP
