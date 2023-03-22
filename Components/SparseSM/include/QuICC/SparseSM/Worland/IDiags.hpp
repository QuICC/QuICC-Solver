/** 
 * @file IDiags.hpp
 * @brief Interface to Worland sparse operator diagonals
 */

#ifndef QUICC_SPARSESM_WORLAND_IDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_IDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/SparseSM/Worland/WorlandKind.hpp"

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
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          * @param q       Truncation q (only consider rows - q equations)
          */
         IDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~IDiags() = default;

         /**
          * @brief Precompute normalization
          *
          * @param maxN Highest radial n
          * @param p    Shift in l
          */
         void precomputeNorm(const int maxN, const int p);
         
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

         /**
          * @brief Zero last n entries of diagonal
          */
         void zeroLast(ACoeff_t& val, const int n) const;

         /**
          * @brief Type of Worland implementation
          */
         Worland::WorlandKind type() const;

         /**
          * @brief Truncation
          */
         const int mQ;

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

         /**
          * @brief Storage for normalization
          */
         std::map<Scalar_t,ACoeff_t> mNorm;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_IDIAGS_HPP
