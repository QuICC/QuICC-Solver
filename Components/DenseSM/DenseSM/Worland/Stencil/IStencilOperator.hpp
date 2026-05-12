/**
 * @file IStencilOperator.hpp
 * @brief Implementation of the generic stencil operator
 */

#ifndef QUICC_DENSESM_WORLAND_STENCIL_ISTENCILOPERATOR_HPP
#define QUICC_DENSESM_WORLAND_STENCIL_ISTENCILOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

/**
 * @brief Implementation of the generic dense stencil operator
 * dense operator
 */
class IStencilOperator : public IWorlandOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of columns
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter: beta = l + dBeta
    * @param l       Harmonic degree l
    */
   IStencilOperator(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l);

   /**
    * @brief Destructor
    */
   virtual ~IStencilOperator() = default;

   /**
    * @brief Get sparse matrix
    */
   SparseMatrix spmat() const;

protected:
   /**
    * @brief Harmonic degree
    */
   const int mL;

private:
};

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_STENCIL_ISTENCILOPERATOR_HPP
