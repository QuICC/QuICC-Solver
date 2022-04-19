/** 
 * @file IDiags.hpp
 * @brief Interface to Worland sparse operator diagonals
 */

#ifndef QUICC_SPARSESM_WORLAND_IDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_IDIAGS_HPP

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
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Interface to Worland sparse operator diagonals
    */ 
   class IDiags
   {
      public:
         /// Typedef for scalar 
         typedef internal::MHDFloat Scalar_t;

         /// Typedef for coefficient array
         typedef internal::ACoeff ACoeff_t;

         /**
          * @brief Constructor
          */
         IDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~IDiags();
         
      protected:
         /**
          * @brief Get alpha
          */
         Scalar_t alpha() const;

         /**
          * @brief Get beta
          */
         Scalar_t beta(const Scalar_t l) const;

         /**
          * @brief Get l
          */
         Scalar_t l() const;

         /**
          * @brief Natural log of norm
          *
          * @param n Row indexes of diagonal
          * @param l l
          */
         ACoeff_t lnorm(const ACoeff_t& n, const Scalar_t l) const;

         /**
          * @brief Unit norm
          *
          * @param n Row indexes of diagonal
          */
         ACoeff_t norm(const ACoeff_t& n) const;

         /**
          * @brief Inverse unit norm
          *
          * @param n Row indexes of diagonal
          */
         ACoeff_t invnorm(const ACoeff_t& n) const;

         /**
          * @brief Normalization for diagonal
          *
          * @param n Row indexes of diagonal
          * @param k Diagonal index
          * @param p Shift in l
          */
         ACoeff_t normalizeDiag(const ACoeff_t& n, const int k, const int p = 0) const;

      private:
         /**
          * @brief Natural log of norm for alpha = -0.5
          *
          * @param n Row indexes of diagonal
          * @param l l
          */
         ACoeff_t lnorm_chebyshev(const ACoeff_t& n, const Scalar_t l) const;

         /**
          * @brief Natural log of norm for alpha = 0.0
          *
          * @param n Row indexes of diagonal
          * @param l l
          */
         ACoeff_t lnorm_legendre(const ACoeff_t& n, const Scalar_t l) const;

         /**
          * @brief Alpha
          */
         Scalar_t mAlpha;

         /**
          * @brief beta = l + dBeta
          */
         Scalar_t mDBeta;

         /**
          * @brief l
          */
         Scalar_t mL;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_IDIAGS_HPP
