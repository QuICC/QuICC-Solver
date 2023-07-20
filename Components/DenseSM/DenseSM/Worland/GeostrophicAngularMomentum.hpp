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
#include "QuICC/Typedefs.hpp"
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
          * @param ugDBeta Geostrophic basis Jacobi beta = l + dBeta
          * @param rows    Number of row
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param q       Truncation q (only consider rows - q equations)
          */
         GeostrophicAngularMomentum(const Scalar_t ugAlpha, const Scalar_t ugDBeta, const int nR, const Scalar_t alpha, const Scalar_t dBeta, const int q = 0);

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
         void buildOpImpl(internal::Matrix& mat, const int rows, const int cols) const final;

      private:
         /**
          * @brief Build Chebyshev type operator
          */
         void buildChebyshevOp(internal::Matrix& mat, const int rows) const;
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_GEOSTROPHICANGULARMOMENTUM_HPP
