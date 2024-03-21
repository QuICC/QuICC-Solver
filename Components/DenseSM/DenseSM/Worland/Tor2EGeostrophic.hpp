/**
 * @file Tor2EGeostrophic.hpp
 * @brief Implementation of the projection operator from the toroidal scalar to
 * the geostrophic basis
 */

#ifndef QUICC_DENSESM_WORLAND_TOR2EGEOSTROPHIC_HPP
#define QUICC_DENSESM_WORLAND_TOR2EGEOSTROPHIC_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IEmbeddedOperator.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the projection operator from the toroidal scalar to
 * the geostrophic basis
 */
class Tor2EGeostrophic : public IEmbeddedOperator
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
    * @param alpha   Geostrophic basis Jacobi alpha
    * @param beta    Geostrophic basis Jacobi beta
    * @param isGenericBasis   Is special geostrophic basis (ie. Li et al)
    * @param alphaB  Jacobi alpha
    * @param betaB   Jacobi beta
    */
   Tor2EGeostrophic(const int nN, const int nL, const int nS, const int nZ,
      const int maxNug, const ArrayI& nli, const int nCpu, const Scalar_t alpha,
      const Scalar_t beta, const bool isGenericBasis, const Scalar_t alphaB,
      const Scalar_t betaB, const bool isTriangular);

   /**
    * @brief Destructor
    */
   virtual ~Tor2EGeostrophic() = default;

protected:
   /**
    * @brief Implementation of build dense matrix operator
    *
    * @param mat operator
    * @param rows rows of matrix
    * @param cols cols of matrix
    */
   void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const final;

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
    * @brief Geostrophic basis alpha
    */
   Scalar_t mUgAlpha;

   /**
    * @brief Geostrophic basis beta
    */
   Scalar_t mUgBeta;

   const bool mIsGenericBasis;

   Scalar_t mAlphaB;
   Scalar_t mBetaB;
   bool mIsTriangular;

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_TOR2EGEOSTROPHIC_HPP
