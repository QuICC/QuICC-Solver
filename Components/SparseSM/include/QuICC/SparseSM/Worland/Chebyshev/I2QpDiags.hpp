/**
 * @file I2QpDiags.hpp
 * @brief Interface to I2Qp diagonals for full sphere Worland I2Qp sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I2QPDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I2QPDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/I2QpDiags.hpp"
#include "QuICC/SparseSM/Worland/Chebyshev/I2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland I2Qp sparse operator
    */
   class I2QpDiags: public QuICC::SparseSM::Worland::I2QpDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param l       Harmonic degree
          * @param q       Truncation q
          */
         I2QpDiags(const Scalar_t alpha, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I2QpDiags() = default;

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

      protected:
         /**
          * @brief Correct diagonal k for q = 1 truncation
          *
          * @param val  Diagonal values to correct
          * @param n    Array of n indexes
          * @param k    Diagonal k
          */
         void correctQ1(ACoeff_t& val, const ACoeff_t& n, const int k) const;

      private:
         /**
          * @brief I2 diagonals for correction
          */
         I2Diags mI2;
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I2QPDIAGS_HPP
