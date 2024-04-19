/**
 * @file GeostrophicEnergy.hpp
 * @brief Implementation of the projection operator from the geostrophic basis to Worland
 */

#ifndef QUICC_DENSESM_WORLAND_GEOSTROPHICENERGY_HPP
#define QUICC_DENSESM_WORLAND_GEOSTROPHICENERGY_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "DenseSM/Worland/IGeostrophicOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the projection operator from the geostrophic basis to Worland
    */
   class GeostrophicEnergy: public IGeostrophicOperator
   {
      public:
         /**
          * @brief Constructor with default Worland basis
          *
          * @param nN      Number of radial modes
          * @param nL      Number of harmonic degrees
          * @param ugAlpha Geostrophic basis Jacobi alpha
          * @param ugBeta  Geostrophic basis Jacobi beta
          * @param isGenericBasis   Use generic Worland basis a geostrophic basis?
          */
         GeostrophicEnergy(const int nN, const int nL, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const bool mIsTriangular);

         /**
          * @brief Destructor
          */
         virtual ~GeostrophicEnergy() = default;

      protected:
         /**
          * @brief Implementation of build dense matrix operator
          *
          * @param mat operator
          * @param rows rows of matrix
          * @param cols cols of matrix
          */
         void buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const final;

         /**
          * @brief Max radial truncation
          */
         const int mNn;

         /**
          * @brief Number of harmonic degrees
          */
         const int mNl;

         /**
          * @brief Uses triangular truncation?
          */
         const bool mIsTriangular;

      private:
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_GEOSTROPHICENERGY_HPP
