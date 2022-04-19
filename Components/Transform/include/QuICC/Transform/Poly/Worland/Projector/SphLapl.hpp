/**
 * @file SphLapl.hpp
 * @brief Implementation of the Worland based spherical laplacian projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_SPHLAPL_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_SPHLAPL_HPP

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
    * @brief Implementation of the Worland based spherical laplacian projector
    */
   class SphLapl: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         SphLapl();

         /**
          * @brief Destructor
          */
         virtual ~SphLapl();

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_SPHLAPL_HPP
