/**
 * @file DivR2FD1R1.hpp
 * @brief Implementation of the spectral operator f/r
 */

#ifndef QUICC_DENSESM_WORLAND_DIVR2FD1R1_HPP
#define QUICC_DENSESM_WORLAND_DIVR2FD1R1_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/ITripleHarmonicOperator.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the spectral operator f/r
 */
class DivR2FD1R1 : public ITripleHarmonicOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lF      harmonic degree of f
    * @param mF      harmonic order of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic order
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    */
   DivR2FD1R1(const int rows, const int cols, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn, std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
      const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~DivR2FD1R1() = default;

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

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_DIVR2FD1R1_HPP
