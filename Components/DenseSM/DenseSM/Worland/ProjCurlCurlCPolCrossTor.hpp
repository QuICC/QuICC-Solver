/**
 * @file ProjCurlCurlCPolCrossTor.hpp
 * @brief Implementation of the r Curl Curl(CurlPolA ^ TorB) projection
 */

#ifndef QUICC_DENSESM_WORLAND_PROJCURLCURLCPOLCROSSTOR_HPP
#define QUICC_DENSESM_WORLAND_PROJCURLCURLCPOLCROSSTOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IProjCrossOperator.hpp"
#include "DenseSM/Worland/ProjCurlCurlTorCrossTor.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the r Curl Curl(CurlPolA ^ TorB) projection
 */
class ProjCurlCurlCPolCrossTor : public ProjCurlCurlTorCrossTor
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param q       Order of quasi-inverse
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lA      harmonic degree of A
    * @param mA      harmonic order of A
    * @param lB      harmonic degree of B
    * @param mB      harmonic order of B
    * @param pPolA   Poloidal radial function
    * @param pTorB   Toroidal radial function
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    */
   ProjCurlCurlCPolCrossTor(const int rows, const int cols, const int q, const int lOut,
      const int mOut, const int lA, const int mA, const int lB, const int mB,
      std::shared_ptr<RadialTorPolFunction> pPolA,
      std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha,
      const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlCurlCPolCrossTor() = default;

protected:
private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_PROJCURLCURLCPOLCROSSTOR_HPP
