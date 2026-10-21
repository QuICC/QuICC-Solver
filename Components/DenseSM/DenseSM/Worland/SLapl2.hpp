/**
 * @file SLapl2.hpp
 * @brief Implementation of the spherical bilaplacian tau operator
 */

#ifndef QUICC_DENSESM_WORLAND_SLAPL2_HPP
#define QUICC_DENSESM_WORLAND_SLAPL2_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/ITauOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the spherical bilaplacian tau operator
 * dense operator
 */
class SLapl2 : public ITauOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of columns
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    * @param l       Harmonic degree l
    */
   SLapl2(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l);

   /**
    * @brief Destructor
    */
   virtual ~SLapl2() = default;

protected:
   /**
    * @brief Implementation of build dense matrix operator
    * @param output operator
    */
   virtual void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const override;
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_SLAPL2_HPP
