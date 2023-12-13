/**
 * @file P.hpp
 * @brief Implementation of the associated Legendre based P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_P_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_P_HPP
// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/IALegendreProjector.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class P;

   /**
    * @brief Implementation of the associated Legendre based P projector
    */
   template<>
   class P<base_t>: public IALegendreProjector
   {
      public:

         /**
          * @brief Constructor
          */
         P() = default;

         /**
          * @brief Destructor
          */
         virtual ~P() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_P_HPP
