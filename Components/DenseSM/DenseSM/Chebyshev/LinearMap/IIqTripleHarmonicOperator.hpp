/**
 * @file IIqTripleHarmonicOperator.hpp
 * @brief Implementation of the generic spectral triple harmonic operator
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IIQTRIPLEEHARMONICOPERATOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IIQTRIPLEEHARMONICOPERATOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ITripleHarmonicOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the generic spectral triple harmonic operator
 */
class IIqTripleHarmonicOperator : public ITripleHarmonicOperator
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
    * @param mOut    Output harmonic degree
    * @param lF      harmonic degree of f
    * @param mF      harmonic degree of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic degree
    * @param pF      Radial function pointer
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   IIqTripleHarmonicOperator(const int rows, const int cols, const int q, const int p, const int lOut,
      const int mOut, const int lF, const int mF, const int lIn, const int mIn,
      std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper) :
       ITripleHarmonicOperator(rows, cols, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower, upper),
       mQ(q) {};

   /**
    * @brief Destructor
    */
   virtual ~IIqTripleHarmonicOperator() = default;

protected:
   /**
    * @brief Order of quasi-inverse
    */
   const int mQ;

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IIQTRIPLEEHARMONICOPERATOR_HPP
