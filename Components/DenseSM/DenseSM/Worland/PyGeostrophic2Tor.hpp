/**
 * @file PyGeostrophic2Tor.hpp
 * @brief Implementation of the projection operator from the geostrophic basis to Worland in Python
 */

#ifndef QUICC_DENSESM_WORLAND_PYGEOSTROPHIC2TOR_HPP
#define QUICC_DENSESM_WORLAND_PYGEOSTROPHIC2TOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the projection operator from the geostrophic basis to Worland in Python
    */
   class PyGeostrophic2Tor: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param ugAlpha Geostrophic basis Jacobi alpha
          * @param ugBeta  Geostrophic basis Jacobi beta
          * @param isGenericBasis   Use generic Worland basis a geostrophic basis?
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         PyGeostrophic2Tor(const int nN, const int nL, const std::vector<int>& nIdx, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const Scalar_t alpha, const Scalar_t dBeta, const bool isTriangular);

         /**
          * @brief Destructor
          */
         virtual ~PyGeostrophic2Tor();

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
          * @brief Max harmonics
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
         /**
          * @brief Build Chebyshev type operator
          */
         void buildChebyshevOp(Internal::Matrix& mat, const int rows, const int cols) const;
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_PYGEOSTROPHIC2TOR_HPP
