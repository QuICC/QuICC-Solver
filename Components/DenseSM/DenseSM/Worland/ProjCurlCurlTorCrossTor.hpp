/**
 * @file ProjCurlCurlTorCrossTor.hpp
 * @brief Implementation of the r Curl Curl(TorA ^ TorB) projection
 */

#ifndef QUICC_DENSESM_WORLAND_PROJCURLCURLTORCROSSTOR_HPP
#define QUICC_DENSESM_WORLAND_PROJCURLCURLTORCROSSTOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IProjCrossOperator.hpp"
#include "DenseSM/Worland/ITripleHarmonicOperator.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the r Curl Curl(TorA ^ TorB) projection
 */
class ProjCurlCurlTorCrossTor : public IProjCrossOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lA      harmonic degree of A
    * @param mA      harmonic order of A
    * @param lB      harmonic degree of B
    * @param mB      harmonic order of B
    * @param pTorA   Toroidal A radial function
    * @param pTorB   Toroidal B radial function
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    */
   ProjCurlCurlTorCrossTor(const int rows, const int cols, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pTorA, std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha, const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlCurlTorCrossTor() = default;

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
    * @brief Radial operator
    */
   std::shared_ptr<ITripleHarmonicOperator>  mpOp;

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_PROJCURLCURLTORCROSSTOR_HPP
