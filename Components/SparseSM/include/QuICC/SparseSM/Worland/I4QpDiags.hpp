/**
 * @file I4QpDiags.hpp
 * @brief Interface to I4Qp diagonals for full sphere Worland I4Qp sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I4QPDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I4QPDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I4Qp sparse operator
    */
   class I4QpDiags: public IDiags
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
         I4QpDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I4QpDiags() = default;

         /**
          * @brief 4. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_4(const ACoeff_t& n) const = 0;

         /**
          * @brief 3. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_3(const ACoeff_t& n) const = 0;

         /**
          * @brief 2. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_2(const ACoeff_t& n) const = 0;

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

         /**
          * @brief 2. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d2(const ACoeff_t& n) const = 0;

         /**
          * @brief 3. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d3(const ACoeff_t& n) const = 0;

         /**
          * @brief 4. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d4(const ACoeff_t& n) const = 0;

         /**
          * @brief 5. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d5(const ACoeff_t& n) const = 0;

      protected:

      private:
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I4QPDIAGS_HPP
