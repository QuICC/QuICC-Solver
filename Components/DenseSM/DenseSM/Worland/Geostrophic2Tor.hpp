/**
 * @file Geostrophic2Tor.hpp
 * @brief Implementation of the projection operator from the geostrophic basis to Worland
 */

#ifndef QUICC_DENSESM_WORLAND_GEOSTROPHIC2TOR_HPP
#define QUICC_DENSESM_WORLAND_GEOSTROPHIC2TOR_HPP

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
   class Geostrophic2Tor: public IGeostrophicOperator
   {
      public:
         /**
          * @brief Constructor with default Worland basis
          *
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugBeta  Geostrophic basis Jacobi beta
          * @param isGenericBasis   Use generic Worland basis a geostrophic basis?
          */
         Geostrophic2Tor(const int nN, const int nL, const std::vector<int>& nIdx, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const bool mIsTriangular);

         /**
          * @brief Constructor
          *
          * @param ugAlpha Geostrophic basis Jacobi alpha
          * @param ugBeta  Geostrophic basis Jacobi beta
          * @param isGenericBasis   Use generic Worland basis a geostrophic basis?
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         Geostrophic2Tor(const int nN, const int nL, const std::vector<int>& nIdx, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const bool isTriangular, const Scalar_t alpha, const Scalar_t dBeta);

         /**
          * @brief Destructor
          */
         virtual ~Geostrophic2Tor() = default;

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
          * @brief radial indexes
          */
         std::vector<int> mNidx;

         /**
          * @brief Uses triangular truncation?
          */
         const bool mIsTriangular;

      private:
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_GEOSTROPHIC2TOR_HPP
