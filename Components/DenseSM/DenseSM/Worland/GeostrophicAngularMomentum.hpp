/**
 * @file GeostrophicAngularMomentum.hpp
 * @brief Implementation of the angular momentum operator for the geostrophic basis
 */

#ifndef QUICC_DENSESM_WORLAND_GEOSTROPHICANGULARMOMENTUM_HPP
#define QUICC_DENSESM_WORLAND_GEOSTROPHICANGULARMOMENTUM_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/Worland/IGeostrophicOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the angula momentum operator for the geostrophic basis
    */
   class GeostrophicAngularMomentum: public IGeostrophicOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugBeta Geostrophic basis Jacobi beta = l + dBeta
          */
         GeostrophicAngularMomentum(const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const int nR);

         /**
          * @brief Constructor
          *
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugBeta Geostrophic basis Jacobi beta = l + dBeta
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         GeostrophicAngularMomentum(const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const int nR, const Scalar_t alpha, const Scalar_t dBeta);

         /**
          * @brief Destructor
          */
         virtual ~GeostrophicAngularMomentum() = default;

      protected:
         /**
          * @brief Implementation of build dense matrix operator
          *
          * @param mat operator
          * @param rows rows of matrix
          * @param cols cols of matrix
          */
         void buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const final;

      private:
         /**
          * @brief Build Worland type independent operator
          */
         void buildGenericOp(Internal::Matrix& mat, const int rows) const;
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_GEOSTROPHICANGULARMOMENTUM_HPP
