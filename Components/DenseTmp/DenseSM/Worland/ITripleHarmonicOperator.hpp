/**
 * @file ITripleHarmonicOperator.hpp
 * @brief Implementation of the generic spectral triple harmonic operator
 */

#ifndef QUICC_DENSESM_WORLAND_ITRIPLEHARMONICOPERATOR_HPP
#define QUICC_DENSESM_WORLAND_ITRIPLEHARMONICOPERATOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the generic spectral triple harmonic operator
 */
class ITripleHarmonicOperator : public IWorlandOperator
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
    * @param pF      F radial function
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    */
   ITripleHarmonicOperator(const int rows, const int cols, const int lOut, const int lF, const int lIn, std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
      const Scalar_t dBeta) : IWorlandOperator(rows, cols, alpha, dBeta), mLout(lOut), mLf(lF), mLin(lIn), mpF(pF){};

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

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_ITRIPLEHARMONICOPERATOR_HPP
