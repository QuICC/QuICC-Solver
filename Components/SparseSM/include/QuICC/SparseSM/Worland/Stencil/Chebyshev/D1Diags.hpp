/**
 * @file D1Diags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland D1 boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_CHEBYSHEV_D1DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_CHEBYSHEV_D1DIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

/// Namespace for Galerkin stencil specific to Worland of Chebyshev kind
namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland D1 boundary condition stencil
    */
   class D1Diags: public QuICC::SparseSM::Worland::Stencil::D1Diags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   jacobi alpha
          * @param l       harmonic degree l
          */
         D1Diags(const Scalar_t alpha, const int l);

         /**
          * @brief Destructor
          */
         virtual ~D1Diags() = default;

         /**
          * @brief 1. subdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d_1(const ACoeff_t& n) const final;

         /**
          * @brief diagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const final;

      protected:

      private:
   };

} // Chebyshev
} // Stencil
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_STENCIL_CHEBYSHEV_D1DIAGS_HPP
