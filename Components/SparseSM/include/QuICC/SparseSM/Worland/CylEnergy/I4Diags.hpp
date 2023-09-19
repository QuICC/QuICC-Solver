/**
 * @file I4Diags.hpp
 * @brief Interface to I4 diagonals for full sphere Worland I4 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CYLENERGY_I4DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CYLENERGY_I4DIAGS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/I4Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace CylEnergy {

   /**
    * @brief Implementation of the full sphere Worland I4 sparse operator
    */
   class I4Diags: public QuICC::SparseSM::Worland::I4Diags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param l       Harmonic degree
          * @param q       Truncation q
          */
         I4Diags(const Scalar_t alpha, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I4Diags() = default;

         /**
          * @brief 4. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_4(const ACoeff_t& n) const;

         /**
          * @brief 3. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_3(const ACoeff_t& n) const;

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

         /**
          * @brief 3. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d3(const ACoeff_t& n) const;

         /**
          * @brief 4. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d4(const ACoeff_t& n) const;

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CYLENERGY_I4DIAGS_HPP
