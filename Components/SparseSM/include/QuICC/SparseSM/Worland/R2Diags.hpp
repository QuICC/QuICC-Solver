/**
 * @file R2Diags.hpp
 * @brief Interface to R2 diagonals for full sphere Worland R2 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_R2DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_R2DIAGS_HPP

// System includes
//

// Project includes
//
#include "Types/Precision.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland R2 sparse operator
    */
   class R2Diags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          * @param q       Truncation q (only consider rows - q equations)
          */
         R2Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~R2Diags() = default;

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

         /**
          * @brief 1. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d1(const ACoeff_t& n) const = 0;

      protected:

      private:
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_R2DIAGS_HPP
