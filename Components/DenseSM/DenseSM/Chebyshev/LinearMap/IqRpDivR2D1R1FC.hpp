/**
 * @file IqRpDivR2D1R1FC.hpp
 * @brief Implementation of the spectral operator I^q r^p 1/r D(r f) (-lapl(*))/r
 *        multiplied by gaunt coefficient
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IQRPDIVR2D1R1FC_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IQRPDIVR2D1R1FC_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ITripleHarmonicOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the spectral operator I^q r^p 1/r D(r f) (-lapl(*))/r
 */
class IqRpDivR2D1R1FC : public ITripleHarmonicOperator
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
    * @param lF      harmonic degree of f
    * @param mF      harmonic order of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic order
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   IqRpDivR2D1R1FC(const int rows, const int cols, const int q,  const int p, const int lOut, const int mOut,
      const int lF, const int mF, const int lIn, const int mIn,
      std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~IqRpDivR2D1R1FC() = default;

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

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IQRPDIVR2D1R1FC_HPP
