/**
 * @file Tor2GridS.hpp
 * @brief Implementation of the projection operator from the toroidal scalar to the geostrophic basis
 */

#ifndef QUICC_DENSESM_WORLAND_TOR2GRIDS_HPP
#define QUICC_DENSESM_WORLAND_TOR2GRIDS_HPP

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
    * @brief Implementation of the projection operator from the toroidal scalar to the geostrophic basis
    */
   class Tor2GridS: public IGeostrophicOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN      Number of radial modes
          * @param nL      Number of harmonic degrees
          * @param nS      Number of cylindrical s modes
          * @param nZ      Number of z grid points
          * @param maxNug  Maximum truncation for geostrophic flow
          * @param nli     Radial truncation nN(l)
          * @param nCpu    Number of CPU in MPI version
          * @param ugAlph  Geostrophic basis Jacobi alpha
          * @param ugDBeta Geostrophic basis Jacobi beta = l + dBeta
          * @param alphaA   Jacobi alpha
          * @param dBetaA   Jacobi beta = l + dBeta
          * @param alphaB   Jacobi alpha
          * @param dBetaB   Jacobi beta = l + dBeta
          */
         Tor2GridS(const int nN, const int nL, const int nS, const int nZ, const int maxNug, const ArrayI& nli, const int nCpu, const Scalar_t ugAlpha, const Scalar_t ugDBeta, const Scalar_t alphaA, const Scalar_t dBetaA, const Scalar_t alphaB, const Scalar_t dBetaB, const bool isTriangular);

         /**
          * @brief Destructor
          */
         virtual ~Tor2GridS() = default;

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
          * @brief Compute quadratue points and weights for z integral
          */
         void computeQuadraturez(Internal::Array& igridz, Internal::Array& iweightz, const int nz) const;

         /**
          * @brief Compute Z integral
          */
         void intgzWorland(int l, int n,
            Internal::Matrix& iintgz, Internal::Array& igrids, Internal::Array& igridz,
            Internal::Array& iweightz) const;

         /**
          * @brief Max radial truncation
          */
         const int mNn;

         /**
          * @brief Number of harmonic degrees
          */
         const int mNl;

         /**
          * @brief Number of cylindrical modes
          */
         const int mNs;

         /**
          * @brief Number of z grid points
          */
         const int mNz;

         /**
          * @brief Number of geostrophic modes
          */
         const int mNnug;

         /**
          * @brief List of radial truncations
          */
         ArrayI mNlist;

         /**
          * @brief Number of CPU
          */
         const int mNcpu;

         /**
          * @brief radial indexes
          */
         std::vector<int> mNidx;

         Scalar_t mAlphaA;
         Scalar_t mDBetaA;
         Scalar_t mAlphaB;
         Scalar_t mDBetaB;
         bool mIsTriangular;

      private:
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_TOR2GRIDS_HPP
