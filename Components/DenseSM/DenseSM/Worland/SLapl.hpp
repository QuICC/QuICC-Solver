/**
 * @file SLapl.hpp
 * @brief Implementation of the spherical laplacian tau operator
 */

#ifndef QUICC_DENSESM_WORLAND_SLAPL_HPP
#define QUICC_DENSESM_WORLAND_SLAPL_HPP

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
 * @brief Implementation of the spherical laplacian tau operator
 * dense operator
 */
class SLapl : public ITauOperator
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
   SLapl(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l);

   /**
    * @brief Destructor
    */
   virtual ~SLapl() = default;

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

#endif // QUICC_DENSESM_WORLAND_SLAPL_HPP
