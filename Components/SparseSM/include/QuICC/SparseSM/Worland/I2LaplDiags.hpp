/**
 * @file I2LaplDiags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland I2Lapl sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I2LAPLDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I2LAPLDIAGS_HPP

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
    * @brief Implementation of the full sphere Worland I2Lapl sparse operator
    */
   class I2LaplDiags: public IDiags
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
         I2LaplDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I2LaplDiags() = default;

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

#endif // QUICC_SPARSESM_WORLAND_I2LAPLDIAGS_HPP
