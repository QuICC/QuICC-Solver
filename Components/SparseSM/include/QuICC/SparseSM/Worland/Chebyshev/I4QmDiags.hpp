/**
 * @file I4QmDiags.hpp
 * @brief Interface to I4Qm diagonals for full sphere Worland I4Qm sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4QMDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4QMDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/I4QmDiags.hpp"
#include "QuICC/SparseSM/Worland/Chebyshev/I4Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland I4Qm sparse operator
    */
   class I4QmDiags: public QuICC::SparseSM::Worland::I4QmDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param l       Harmonic degree
          * @param q       Truncation q
          */
         I4QmDiags(const Scalar_t alpha, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I4QmDiags() = default;

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
          * @brief 2. superdiagonal
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

         /**
          * @brief 5. superdiagonal
          *
          * @param n Array of n indexes
          */
         ACoeff_t d5(const ACoeff_t& n) const final;

      protected:
         /**
          * @brief Correct diagonal k for q = 2 truncation
          *
          * @param val  Diagonal values to correct
          * @param n    Array of n indexes
          * @param k    Diagonal k
          */
         void correctQ2(ACoeff_t& val, const ACoeff_t& n, const int k) const;

      private:
         /**
          * @brief I4 diagonals for correction
          */
         I4Diags mI4;
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4QMDIAGS_HPP
