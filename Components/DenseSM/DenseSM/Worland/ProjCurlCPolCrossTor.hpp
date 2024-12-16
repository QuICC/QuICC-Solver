/**
 * @file ProjCurlCPolCrossTor.hpp
 * @brief Implementation of the r Curl(CurlPolA ^ TorB) projection
 */

#ifndef QUICC_DENSESM_WORLAND_PROJCURLCPOLCROSSTOR_HPP
#define QUICC_DENSESM_WORLAND_PROJCURLCPOLCROSSTOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IProjCrossOperator.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "DenseSM/Worland/ProjCurlTorCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the r Curl(curlPolA ^ TorB) projection
 */
class ProjCurlCPolCrossTor : public ProjCurlTorCrossTor
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
    * @param pPolA   Poloidal A radial function
    * @param pTorB   Toroidal B radial function
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    */
   ProjCurlCPolCrossTor(const int rows, const int cols, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha, const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlCPolCrossTor() = default;

protected:

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_PROJCURLCPOLCROSSTOR_HPP
