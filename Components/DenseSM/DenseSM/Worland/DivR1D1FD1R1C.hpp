/**
 * @file DivR1D1FD1R1C.hpp
 * @brief Implementation of the spectral operator 1/r D(f) D(r lapl *)
 */

#ifndef QUICC_DENSESM_WORLAND_DIVR1D1FD1R1C_HPP
#define QUICC_DENSESM_WORLAND_DIVR1D1FD1R1C_HPP

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
 * @brief Implementation of the spectral operator 1/r D(f (-lapl(*)))
 */
class DivR1D1FD1R1C : public ITripleHarmonicOperator
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
   DivR1D1FD1R1C(const int rows, const int cols, const int lOut, const int mOut,
      const int lF, const int mF, const int lIn, const int mIn,
      std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
      const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~DivR1D1FD1R1C() = default;

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

#endif // QUICC_DENSESM_WORLAND_DIVR1D1FD1R1C_HPP
