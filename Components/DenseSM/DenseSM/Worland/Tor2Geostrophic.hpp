/**
 * @file Tor2Geostrophic.hpp
 * @brief Implementation of the projection operator from the toroidal scalar to
 * the geostrophic basis
 */

#ifndef QUICC_DENSESM_WORLAND_TOR2GEOSTROPHIC_HPP
#define QUICC_DENSESM_WORLAND_TOR2GEOSTROPHIC_HPP

// System includes
//

// Project includes
//
#include "DenseSM/IMatrixSMOperator.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the projection operator from the toroidal scalar to
 * the geostrophic basis
 */
class Tor2Geostrophic : public IMatrixSMOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param nN      Number of radial modes
    * @param nL      Number of harmonic degrees
    * @param nCpu    Number of CPU in MPI version
    * @param alpha   Geostrophic basis Jacobi alpha
    * @param beta    Geostrophic basis Jacobi beta
    * @param isGenericBasis   Is generic Worland basis (ie. Li et al)
    * @param wAlpha  Worland Jacobi alpha
    * @param wDBeta  Worland Jacobi beta = l + dBeta
    */
   Tor2Geostrophic(const int nN, const int nL, const int nCpu, const Scalar_t alpha,
      const Scalar_t beta, const bool isGenericBasis, const bool isTriangular,
      const Scalar_t wAlpha, const Scalar_t wDBeta);

   /**
    * @brief Constructor
    *
    * @param nN      Number of radial modes
    * @param nL      Number of harmonic degrees
    * @param nCpu    Number of CPU in MPI version
    * @param alpha   Geostrophic basis Jacobi alpha
    * @param beta    Geostrophic basis Jacobi beta
    * @param isGenericBasis   Is generic Worland basis (ie. Li et al)
    */
   Tor2Geostrophic(const int nN, const int nL, const int nCpu, const Scalar_t alpha,
      const Scalar_t beta, const bool isGenericBasis, const bool isTriangular);

   /**
    * @brief Destructor
    */
   virtual ~Tor2Geostrophic() = default;

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

   /**
    * @brief Use generic Worland basis as Geostrophic basis?
    */
   const bool mIsGenericBasis;

   bool mIsTriangular;

   /**
    * @brief Worland Jacobi alpha
    */
   Scalar_t mAlpha;

   /**
    * @brief Worland Jacobi beta = l + dBeta
    */
   Scalar_t mDBeta;

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_TOR2GEOSTROPHIC_HPP
