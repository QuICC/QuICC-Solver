/**
 * @file ProjCurlPolCrossPol.hpp
 * @brief Implementation of the r Curl(PolA ^ PolB) projection
 */

#ifndef QUICC_DENSESM_WORLAND_PROJCURLPOLCROSSPOL_HPP
#define QUICC_DENSESM_WORLAND_PROJCURLPOLCROSSPOL_HPP

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
 * @brief Implementation of the r Curl(PolA ^ PolB) projection
 */
class ProjCurlPolCrossPol : public IProjCrossOperator
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
    * @param mB      harmonic order of B
    * @param pPolA   Poloidal A radial function
    * @param pPolB   Poloidal B radial function
    * @param alpha   Jacobi alpha parametery
    * @param dBeta   Jacobi dBeta parameter
    */
   ProjCurlPolCrossPol(const int rows, const int cols, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t alpha, const Scalar_t dBeta);

   /**
    * @brief Destructor
    */
   virtual ~ProjCurlPolCrossPol() = default;

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
   std::shared_ptr<ITripleHarmonicOperator>  mpOpA;

   /**
    * @brief Radial operator B
    */
   std::shared_ptr<ITripleHarmonicOperator>  mpOpB;

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_PROJCURLPOLCROSSPOL_HPP
