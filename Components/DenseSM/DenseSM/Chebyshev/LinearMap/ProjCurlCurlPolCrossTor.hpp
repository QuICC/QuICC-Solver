/**
 * @file ProjCurlCurlPolCrossTor.hpp
 * @brief Implementation of the r Curl Curl(PolA ^ TorB) projection
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLPOLCROSSTOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLPOLCROSSTOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/IProjCrossOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/ITripleHarmonicOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the r Curl Curl(PolA ^ TorB) projection
 */
class ProjCurlCurlPolCrossTor : public IProjCrossOperator
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
   ProjCurlCurlPolCrossTor(const int rows, const int cols, const int p,
      const int lOut, const int mOut, const int lA, const int mA, const int lB,
      const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
      std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlCurlPolCrossTor() = default;

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
    * @brief Radial operator A
    */
   std::shared_ptr<ITripleHarmonicOperator> mpOpA;

   /**
    * @brief Radial operator B
    */
   std::shared_ptr<ITripleHarmonicOperator> mpOpB;

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLPOLCROSSTOR_HPP
