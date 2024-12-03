/**
 * @file ITripleHarmonicOperator.hpp
 * @brief Implementation of the generic spectral triple harmonic operator
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ITRIPLEHARMONICOPERATOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ITRIPLEHARMONICOPERATOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the generic spectral triple harmonic operator
 */
class ITripleHarmonicOperator : public ILinearMapOperator
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
   ITripleHarmonicOperator(const int rows, const int cols, const int lOut, const int lF, const int lIn, std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper) : ILinearMapOperator(rows, cols, lower, upper), mLout(lOut), mLf(lF), mLin(lIn), mpF(pF){};

   /**
    * @brief Destructor
    */
   virtual ~ITripleHarmonicOperator() = default;

protected:
   /**
    * @brief Harmonic degree of output
    */
   int mLout;

   /**
    * @brief Harmonic degree of f
    */
   int mLf;

   /**
    * @brief Harmonic degree of input
    */
   int mLin;

   /**
    * @brief Functor for f
    */
   std::shared_ptr<RadialTorPolFunction> mpF;

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ITRIPLEHARMONICOPERATOR_HPP
