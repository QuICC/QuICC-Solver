/**
 * @file Llm1DivS1.hpp
 * @brief Implementation of the associated Legendre based 1/sin [l(l+1)-1] P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LLM1DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LLM1DIVS1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class Llm1DivS1;

   /**
    * @brief Implementation of the associated Legendre based 1/sin [l(l+1)-1] P projector
    */
   template <>
   class Llm1DivS1<base_t>: public DivS1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         Llm1DivS1() = default;

         /**
          * @brief Destructor
          */
         virtual ~Llm1DivS1() = default;

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
         mutable Array mLlm1;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_BASE_LLDIVS1_HPP
