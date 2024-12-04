/**
 * @file R2DivR2FD1R1.hpp
 * @brief Implementation of the spectral operator r^2 1/r^2 f D(r *)
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_R2DIVR2FD1R1_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_R2DIVR2FD1R1_HPP

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
 * @brief Implementation of the spectral operator r^2 1/r^2 f D(r *)
 */
class R2DivR2FD1R1 : public ITripleHarmonicOperator
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
   R2DivR2FD1R1(const int rows, const int cols, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn, std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~R2DivR2FD1R1() = default;

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

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_R2DIVR2FD1R1_HPP
