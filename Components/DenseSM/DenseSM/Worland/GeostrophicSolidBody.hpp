/**
 * @file GeostrophicSolidBody.hpp
 * @brief Implementation of the solid body with unit angular momentum operator for the geostrophic basis
 */

#ifndef QUICC_DENSESM_WORLAND_GEOSTROPHICSOLIDBODY_HPP
#define QUICC_DENSESM_WORLAND_GEOSTROPHICSOLIDBODY_HPP

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
    * @brief Implementation of the solid body with unit angular momentum operator for the geostrophic basis
    */
   class GeostrophicSolidBody: public IGeostrophicOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN      Number of radial modes
          * @param ugAlpha Geostrophic basis Jacobi alpha
          * @param ugBeta  Geostrophic basis Jacobi beta = l + dBeta
          * @param isGenericBasis   Use generic Worland basis a geostrophic basis?
          */
         GeostrophicSolidBody(const int nN, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis);

         /**
          * @brief Destructor
          */
         virtual ~GeostrophicSolidBody() = default;

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

#endif // QUICC_DENSESM_WORLAND_GEOSTROPHICSOLIDBODY_HPP
