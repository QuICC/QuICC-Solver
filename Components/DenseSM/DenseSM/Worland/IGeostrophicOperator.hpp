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
#include "DenseSM/Worland/IEmbeddedOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the base for geostrophic basis operator
    */
   class IGeostrophicOperator: public IEmbeddedOperator
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
