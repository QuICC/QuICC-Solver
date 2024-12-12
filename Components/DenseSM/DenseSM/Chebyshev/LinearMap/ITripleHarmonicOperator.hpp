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
    * @param mOut    Output harmonic degree
    * @param lF      harmonic degree of f
    * @param mF      harmonic degree of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic degree
    * @param pF      Radial function pointer
    * @param lower   Lower boundar
    * @param upper   Upper boundar
    */
   ITripleHarmonicOperator(const int rows, const int cols, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn, std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper) : ILinearMapOperator(rows, cols, lower, upper), mLout(lOut), mMout(mOut), mLf(lF), mMf(mF), mLin(lIn), mMin(mIn), mpF(pF){};

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
    * @brief Harmonic order of output
    */
   int mMout;

   /**
    * @brief Harmonic degree of f
    */
   int mLf;

   /**
    * @brief Harmonic order of f
    */
   int mMf;

   /**
    * @brief Harmonic degree of input
    */
   int mLin;

   /**
    * @brief Harmonic order of input
    */
   int mMin;

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
