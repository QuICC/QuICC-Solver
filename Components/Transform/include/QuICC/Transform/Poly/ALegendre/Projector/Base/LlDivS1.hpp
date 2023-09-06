/**
 * @file LlDivS1.hpp
 * @brief Implementation of the associated Legendre based 1/sin l(l+1) P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LLDIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LLDIVS1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class LlDivS1;

   /**
    * @brief Implementation of the associated Legendre based 1/sin l(l+1) P projector
    */
   template <>
   class LlDivS1<base_t>: public DivS1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         LlDivS1() = default;

         /**
          * @brief Destructor
          */
         virtual ~LlDivS1() = default;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LLDIVS1_HPP
