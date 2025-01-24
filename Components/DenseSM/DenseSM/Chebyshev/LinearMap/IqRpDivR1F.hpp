/**
 * @file IqRpDivR1F.hpp
 * @brief Implementation of the spectral operator I^q r^p f/r
 *        multiplied by elsasser coefficient
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IQRPDIVR1F_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IQRPDIVR1F_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/IIqTripleHarmonicOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the spectral operator I^q r^p f/r
 */
class IqRpDivR1F : public IIqTripleHarmonicOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param q       Order of quasi-inverse
    * @param p       Power for radial prefactor
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lF      harmonic degree of f
    * @param mF      harmonic order of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic order
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   IqRpDivR1F(const int rows, const int cols, const int q, const int p, const int lOut, const int mOut,
      const int lF, const int mF, const int lIn, const int mIn,
      std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~IqRpDivR1F() = default;

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

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IQRPDIVR1F_HPP
