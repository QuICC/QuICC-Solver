/**
 * @file WorlandKind.hpp
 * @brief Existing Worland basis kinds
 */

#ifndef QUICC_DENSESM_WORLAND_WORLANDKIND_HPP
#define QUICC_DENSESM_WORLAND_WORLANDKIND_HPP

// Project includes
//

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /// Different kinds of Worland polynomials
   enum class WorlandKind {
      CHEBYSHEV,
      LEGENDRE,
      CYLENERGY,
      SPHENERGY,
   };

}
}
}

#endif // QUICC_DENSESM_WORLAND_WORLANDKIND_HPP
