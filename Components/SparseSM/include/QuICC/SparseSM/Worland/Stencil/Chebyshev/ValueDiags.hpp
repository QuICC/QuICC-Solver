/**
 * @file ValueDiags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland Value boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_CHEBYSHEV_VALUEDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_CHEBYSHEV_VALUEDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland Value boundary condition stencil
    */
   class ValueDiags: public QuICC::SparseSM::Worland::Stencil::ValueDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   jacobi alpha
          * @param l       harmonic degree l
          */
         ValueDiags(const Scalar_t alpha, const int l);

         /**
          * @brief Destructor
          */
         virtual ~ValueDiags() = default;

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

#endif // QUICC_SPARSESM_WORLAND_STENCIL_CHEBYSHEV_VALUEDIAGS_HPP
