/**
 * @file ProjCurlPolCrossTor.hpp
 * @brief Implementation of the r Curl(PolA ^ TorB) projection
 */

#ifndef QUICC_DENSESM_WORLAND_PROJCURLPOLCROSSTOR_HPP
#define QUICC_DENSESM_WORLAND_PROJCURLPOLCROSSTOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IProjCrossOperator.hpp"
#include "DenseSM/Worland/ITripleHarmonicOperator.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the r Curl(PolA ^ TorB) projection
 */
class ProjCurlPolCrossTor : public IProjCrossOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lA      harmonic degree of A
    * @param mA      harmonic order of A
    * @param lB      harmonic degree of B
    * @param mB      harmonic ordera of B
    * @param pPolA   Poloidal A radial function
    * @param pTorB   Toroidal B radial function
    * @param alpha   Jacobi alpha parametery
    * @param dBeta   Jacobi dBeta parameter
    */
   ProjCurlPolCrossTor(const int rows, const int cols, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha, const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlPolCrossTor() = default;

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
    * @brief Radial operator
    */
   std::shared_ptr<ITripleHarmonicOperator>  mpOp;

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_PROJCURLPOLCROSSTOR_HPP
