/**
 * @file WorlandKind.hpp
 * @brief Existing Worland basis kinds
 */

#ifndef QUICC_SPARSESM_WORLAND_WORLANDKIND_HPP
#define QUICC_SPARSESM_WORLAND_WORLANDKIND_HPP

// Project includes
//

namespace QuICC {

namespace SparseSM {

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

#endif // QUICC_SPARSESM_WORLAND_WORLANDKIND_HPP
