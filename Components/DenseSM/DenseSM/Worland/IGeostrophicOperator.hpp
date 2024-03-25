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
#include "DenseSM/Worland/IWorlandOperator.hpp"

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
          * @brief Constructor with default Worland basis
          *
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugBeta  Geostrophic basis Jacobi beta
          * @param isGenericBasis   Use generic Worland as geostrophic basis
          * @param rows    Number of row
          * @param cols    Number of cols
          */
         IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugBeta, bool isGenericBasis, const int rows, const int cols);

         /**
          * @brief Constructor
          *
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugBeta  Geostrophic basis Jacobi beta
          * @param isGenericBasis   Use generic Worland as geostrophic basis
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param alpha   Worland Jacobi alpha
          * @param dBeta   Worland Jacobi beta = l + dBeta
          */
         IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugBeta, bool isGenericBasis, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta);

         /**
          * @brief Destructor
          */
         virtual ~IGeostrophicOperator() = default;

      protected:
         /**
          * @brief Geostrophic alpha
          */
         const Scalar_t mcUgAlpha;

         /**
          * @brief Geostrophic beta
          */
         const Scalar_t mcUgBeta;

         /**
          * @brief Use generic Worland basis for Geostrophic basis?
          */
         const bool mcIsGenericBasis;

      private:
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_IGEOSTROPHICOPERATOR_HPP
