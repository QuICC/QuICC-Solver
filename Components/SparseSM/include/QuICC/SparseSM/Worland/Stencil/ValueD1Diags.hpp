/**
 * @file ValueD1Diags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland ValueD1 boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_VALUED1DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_VALUED1DIAGS_HPP

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
    * @brief Implementation of the full sphere Worland ValueD1 boundary condition stencil
    */
   class ValueD1Diags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         ValueD1Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~ValueD1Diags() = default;

         /**
          * @brief 2. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_2(const ACoeff_t& n) const = 0;

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

#endif // QUICC_SPARSESM_WORLAND_STENCIL_VALUED1DIAGS_HPP
