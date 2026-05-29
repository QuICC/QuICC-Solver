/**
 * @file ProjPolViscD2.hpp
 * @brief Implementation of the poloidal projection of the term -D2 u /rho
 */

#ifndef QUICC_DENSESM_WORLAND_PROJPOLVISCD2_HPP
#define QUICC_DENSESM_WORLAND_PROJPOLVISCD2_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/ITripleHarmonicOperator.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the spectral operator -D2 u /rho
 */
class ProjPolViscD2 : public ITripleHarmonicOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param nNr     Number of rows
    * @param nNc     Number of cols
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lF      harmonic degree of f
    * @param mF      harmonic order of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic order
    * @param pF      F radial function
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    */
   ProjPolViscD2(const int nNr, const int nNc, const int lOut, const int mOut,
      const int lF, const int mF, const int lIn, const int mIn,
      std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
      const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~ProjPolViscD2() = default;

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

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_PROJPOLVISCD2_HPP
