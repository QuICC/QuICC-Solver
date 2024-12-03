/**
 * @file DivR1F.hpp
 * @brief Implementation of the spectral operator f/r
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_DIVR1F_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_DIVR1F_HPP

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
 * @brief Implementation of the spectral operator f/r
 */
class DivR1F : public ITripleHarmonicOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param lOut    Output harmonic degree
    * @param lF      harmonic degree of f
    * @param lIn     Input harmonic degree
    * @param lower   Lower boundar
    * @param upper   Upper boundar
    */
   DivR1F(const int rows, const int cols, const int lOut, const int lF, const int lIn, std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~DivR1F() = default;

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

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_DIVR1F_HPP
