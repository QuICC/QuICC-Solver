/**
 * @file ProjCurlCurlCPolCrossTor.hpp
 * @brief Implementation of the r Curl Curl(CurlPolA ^ TorB) projection
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLCPOLCROSSTOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLCPOLCROSSTOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/IProjCrossOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCurlTorCrossTor.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

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
    * @param p       Power of radial prefactor
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lA      harmonic degree of A
    * @param mA      harmonic order of A
    * @param lB      harmonic degree of B
    * @param mB      harmonic order of B
    * @param pPolA   Poloidal radial function
    * @param pTorB   Toroidal radial function
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   ProjCurlCurlCPolCrossTor(const int rows, const int cols, const int q, const int p,
      const int lOut, const int mOut, const int lA, const int mA, const int lB,
      const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
      std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlCurlCPolCrossTor() = default;

protected:
private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLCPOLCROSSTOR_HPP
