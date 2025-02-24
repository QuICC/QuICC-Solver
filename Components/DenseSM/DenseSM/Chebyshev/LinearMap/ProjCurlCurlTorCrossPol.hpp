/**
 * @file ProjCurlCurlTorCrossPol.hpp
 * @brief Implementation of the r Curl Curl(TorA ^ PolB) projection
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLTORCROSSPOL_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLTORCROSSPOL_HPP

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
 * @brief Implementation of the r Curl Curl(TorA ^ PolB) projection
 */
class ProjCurlCurlTorCrossPol : public IProjCrossOperator
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
    * @param pTorA   Toroidal radial function
    * @param pPolB   Toroidal radial function
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   ProjCurlCurlTorCrossPol(const int rows, const int cols, const int q, const int p,
      const int lOut, const int mOut, const int lA, const int mA, const int lB,
      const int mB, std::shared_ptr<RadialTorPolFunction> pTorA,
      std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlCurlTorCrossPol() = default;

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

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_PROJCURLCURLTORCROSSPOL_HPP
