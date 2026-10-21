/**
 * @file ITauOperator.hpp
 * @brief Implementation of the generic tau operator
 */

#ifndef QUICC_DENSESM_WORLAND_ITAUOPERATOR_HPP
#define QUICC_DENSESM_WORLAND_ITAUOPERATOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the generic dense stencil operator
 * dense operator
 */
class ITauOperator : public IWorlandOperator
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
   ITauOperator(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l);

   /**
    * @brief Destructor
    */
   virtual ~ITauOperator() = default;

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

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_ITAUOPERATOR_HPP
