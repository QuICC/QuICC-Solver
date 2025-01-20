/**
 * @file ProjCurlCPolCrossTor.hpp
 * @brief Implementation of the r Curl(CurlPolA ^ TorB) projection
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCPOLCROSSTOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCPOLCROSSTOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/IProjCrossOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "DenseSM/Chebyshev/LinearMap/ProjCurlTorCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

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
    * @param p       Power of radial prefactor
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lA      harmonic degree of A
    * @param mA      harmonic order of A
    * @param lB      harmonic degree of B
    * @param mB      harmonic order of B
    * @param pPolA   Poloidal A radial function
    * @param pTorB   Toroidal B radial function
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   ProjCurlCPolCrossTor(const int rows, const int cols, const int p, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower, const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlCPolCrossTor() = default;

protected:

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCPOLCROSSTOR_HPP
