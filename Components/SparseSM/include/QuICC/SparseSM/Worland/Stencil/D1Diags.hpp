/**
 * @file D1Diags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland D1 boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_D1DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_D1DIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

/// Namespace for Galerkin stencil matrices for Worland polynomials
namespace Stencil {

   /**
    * @brief Implementation of the full sphere Worland D1 boundary condition stencil
    */
   class D1Diags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         D1Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~D1Diags() = default;

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

#endif // QUICC_SPARSESM_WORLAND_STENCIL_D1DIAGS_HPP
