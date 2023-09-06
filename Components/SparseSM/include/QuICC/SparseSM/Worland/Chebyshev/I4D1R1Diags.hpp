/**
 * @file I4D1R1Diags.hpp
 * @brief Interface to I4D1R1 diagonals for full sphere Worland I4D1R1 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4D1R1DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4D1R1DIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/I4D1R1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland I4D1R1 sparse operator
    */
   class I4D1R1Diags: public QuICC::SparseSM::Worland::I4D1R1Diags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param l       Harmonic degree
          * @param q       Truncation q
          */
         I4D1R1Diags(const Scalar_t alpha, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I4D1R1Diags() = default;

         /**
          * @brief 3. subdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d_3(const ACoeff_t& n) const final;

         /**
          * @brief 2. subdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d_2(const ACoeff_t& n) const final;

         /**
          * @brief 1. subdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d_1(const ACoeff_t& n) const final;

         /**
          * @brief Main diagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const final;

         /**
          * @brief 1. superdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d1(const ACoeff_t& n) const final;

         /**
          * @brief 2. superdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d2(const ACoeff_t& n) const final;

         /**
          * @brief 3. superdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d3(const ACoeff_t& n) const final;

         /**
          * @brief 4. superdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d4(const ACoeff_t& n) const final;

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4D1R1DIAGS_HPP
