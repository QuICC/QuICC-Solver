/**
 * @file R2Diags.hpp
 * @brief Interface to R2 diagonals for full sphere Worland R2 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_R2DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_R2DIAGS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/R2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland R2 sparse operator
    */
   class R2Diags: public QuICC::SparseSM::Worland::R2Diags
   {
      public:
         /**
          * @brief Constructor
          */
         R2Diags(const Scalar_t alpha, const int l);

         /**
          * @brief Destructor
          */
         virtual ~R2Diags() = default;

         /**
          * @brief 1. subdiagonal
          */
         ACoeff_t d_1(const ACoeff_t& n) const final;

         /**
          * @brief Main diagonal
          */
         ACoeff_t d0(const ACoeff_t& n) const final;

         /**
          * @brief 1. superdiagonal
          */
         ACoeff_t d1(const ACoeff_t& n) const final;

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_R2DIAGS_HPP
