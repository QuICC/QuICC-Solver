/**
 * @file IGeostrophicOperator.hpp
 * @brief Implementation of the base for geostrophic basis operator
 */

#ifndef QUICC_DENSESM_WORLAND_IGEOSTROPHICOPERATOR_HPP
#define QUICC_DENSESM_WORLAND_IGEOSTROPHICOPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "DenseSM/IWorlandOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the base for geostrophic basis operator
    */
   class IGeostrophicOperator: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugDBeta Geostrophic basis Jacobi beta = l + dBeta
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param q       Truncation q (only consider rows - q equations)
          */
         IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugDBeta, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int q = 0);

         /**
          * @brief Destructor
          */
         virtual ~IGeostrophicOperator() = default;

      protected:
         /**
          * @brief Is Ug basis?
          */
         bool isUgBasis() const;

         /**
          * @brief alpha and beta parameter of Ug basis?
          */
         bool isUgBasis(const Scalar_t a, const Scalar_t b) const;

         /**
          * @brief Akj coefficient used for geostrophic to toroidal projection (Jiawen's thesis: (C.3))
          */
         Internal::MHDFloat Akj(const int k, const int j) const;

         /**
          * @brief Bjn coefficient used for geostrophic to toroidal projection (Jiawen's thesis: (3.28) but using Cn in place of Cnab)
          */
         Internal::MHDFloat Bjn(const int j, const int n) const;

         /**
          * @brief Bjnab coefficient used for geostrophic to toroidal projection (Jiawen's thesis: (3.28))
          */
         Internal::MHDFloat Bjnab(const int j, const int n, const Internal::MHDFloat a, const Internal::MHDFloat b) const;

         /**
          * @Brief Normalization coefficient for $\Lambda_n(s)$ (Jiawen's thesis (3.20))
          */
         Internal::MHDFloat Cnab(const int n, const Scalar_t a, const Scalar_t b) const;

         /**
          * @Brief Normalization coefficient for $\tilde{\Lambda}_n(s)$ (Jiawen's thesis (3.23))
          */
         Internal::MHDFloat Cn(const int n) const;

         /**
          * @brief Geostrophic alpha
          */
         Scalar_t mcUgAlpha;

         /**
          * @brief Geostrophic dBeta: beta = l + dbeta
          */
         Scalar_t mcUgDBeta;

      private:
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_IGEOSTROPHICOPERATOR_HPP
