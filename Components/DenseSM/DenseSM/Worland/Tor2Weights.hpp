/**
 * @file Tor2Weights.hpp
 * @brief Implementation of the projection operator from the toroidal scalar to
 * the geostrophic basis
 */

#ifndef QUICC_DENSESM_WORLAND_TOR2WEIGHTS_HPP
#define QUICC_DENSESM_WORLAND_TOR2WEIGHTS_HPP

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
class Tor2Weights : public IMatrixSMOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param nN      Number of radial modes
    * @param nL      Number of harmonic degrees
    * @param alpha   Jacobi alpha
    * @param beta    Jacobi beta
    * @param isTriangular  Uses triangular truncation?
    */
   Tor2Weights(const int nN, const int nL, const Scalar_t alpha,
      const Scalar_t beta, const bool isTriangular);

   /**
    * @brief Destructor
    */
   virtual ~Tor2Weights() = default;

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
    * @brief Jacobi alpha
    */
   Scalar_t mAlpha;

   /**
    * @brief Jacobi beta
    */
   Scalar_t mBeta;

   /**
    * @brief Uses triangular truncation?
    */
   bool mIsTriangular;

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_TOR2WEIGHTS_HPP
