/**
 * @file I2Diags.hpp
 * @brief Interface to I2 diagonals for full sphere Worland I2 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CYLENERGY_I2DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CYLENERGY_I2DIAGS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/I2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

/// Namespace for sparse operator for Worland polynomials of CylEnergy kind (orthogonal with cylindrical energy weight)
namespace CylEnergy {

   /**
    * @brief Implementation of the full sphere Worland I2 sparse operator
    */
   class I2Diags: public QuICC::SparseSM::Worland::I2Diags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param l       Harmonic degree
          * @param q       Truncation q
          */
         I2Diags(const Scalar_t alpha, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I2Diags() = default;

         /**
          * @brief 2. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_2(const ACoeff_t& n) const;

         /**
          * @brief 1. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_1(const ACoeff_t& n) const;

         /**
          * @brief Main diagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const;

         /**
          * @brief 1. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d1(const ACoeff_t& n) const;

         /**
          * @brief 2. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d2(const ACoeff_t& n) const;

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CYLENERGY_I2DIAGS_HPP
