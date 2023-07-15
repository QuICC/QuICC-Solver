/**
 * @file Ll.hpp
 * @brief Implementation of the associated Legendre based l(l+1) P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LL_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class Ll;

   /**
    * @brief Implementation of the associated Legendre based l(l+1) P projector
    */
   template <>
   class Ll<base_t>: public P<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         Ll() = default;

         /**
          * @brief Destructor
          */
         virtual ~Ll() = default;

      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;

         /**
          * @brief l(l+1) scaling factors
          */
         virtual void initSpecial() const;

         /**
          * @brief Storage for l dependent factors
          */
         mutable Array mLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LL_HPP
