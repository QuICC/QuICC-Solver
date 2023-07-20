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
#include "QuICC/Typedefs.hpp"
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
          * @brief Constructor
          *
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugDBeta Geostrophic basis Jacobi beta = l + dBeta
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param q       Truncation q (only consider rows - q equations)
          */
         Geostrophic2Tor(const int nN, const int maxnl, const int nR, const int maxNug, const ArrayI& nli, const std::vector<int>& nIdx, const Scalar_t ugAlpha, const Scalar_t ugDBeta, const Scalar_t alpha, const Scalar_t dBeta, const int q = 0);

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
         void buildOpImpl(internal::Matrix& mat, const int rows, const int cols) const final;

         /**
          * @brief Max radial truncation
          */
         const int mNn;

         /**
          * @brief Max harmonics
          */
         const int mMaxnl;

         /**
          * @brief radial grid size
          */
         const int mNr;

         /**
          * @brief Max geostrophic truncation
          */
         const int mMaxNug;

         /**
          * @brief List of radial truncations
          */
         ArrayI mNlist;

         /**
          * @brief radial indexes
          */
         std::vector<int> mNidx;

      private:
         /**
          * @brief Build Chebyshev type operator
          */
         void buildChebyshevOp(internal::Matrix& mat, const int rows, const int cols) const;
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_GEOSTROPHIC2TOR_HPP
