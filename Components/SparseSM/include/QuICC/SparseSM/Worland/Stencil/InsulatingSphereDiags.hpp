/**
 * @file InsulatingSphereDiags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland insulating sphere boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_INSULATINGSPHEREDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_INSULATINGSPHEREDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   /**
    * @brief Implementation of the full sphere Worland insulating sphere boundary condition stencil
    */
   class InsulatingSphereDiags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         InsulatingSphereDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~InsulatingSphereDiags() = default;

         /**
          * @brief 1. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_1(const ACoeff_t& n) const = 0;

         /**
          * @brief Main diagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const = 0;

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_STENCIL_INSULATINGSPHEREDIAGS_HPP
