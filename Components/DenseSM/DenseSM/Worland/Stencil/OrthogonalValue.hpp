/**
 * @file OrthogonalValue.hpp
 * @brief Implementation of the orthogonal Galerkin stencil for Value boundary condition
 *        based on: https://homepages.see.leeds.ac.uk/~earpwl/Galerkin/Galerkin.html
 *        Livermore, P., 2010. Galerkin orthogonal polynomials, J. Comp. Phys., 229(6), 2046–2060.
 */

#ifndef QUICC_DENSESM_WORLAND_STENCIL_ORTHOGONALVALUE_HPP
#define QUICC_DENSESM_WORLAND_STENCIL_ORTHOGONALVALUE_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "DenseSM/Worland/Stencil/IStencilOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

/**
 * @brief Implementation of the orthogonal Galerkin stencil for Value boundary condition
 * dense operator
 */
class OrthogonalValue : public IStencilOperator
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
   OrthogonalValue(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l);

   /**
    * @brief Destructor
    */
   virtual ~OrthogonalValue() = default;

protected:
   /**
    * @brief Implementation of build dense matrix operator
    * @param output operator
    */
   virtual void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const override;

private:
   /**
    * @brief c1 expansion coefficient
    */
   Scalar_t c1(const int n, const int l) const;

   /**
    * @brief c2 expansion coefficient
    */
   Scalar_t c2(const int n, const int l) const;

   /**
    * @brief c3 expansion coefficient
    */
   Scalar_t c3(const int n, const int l) const;

   /**
    * @brief rescaled c1 expansion coefficient
    */
   Scalar_t d1(const int n, const int l) const;

   /**
    * @brief rescaled c2 expansion coefficient
    */
   Scalar_t d2(const int n, const int l) const;

   /**
    * @brief rescaled c3 expansion coefficient
    */
   Scalar_t d3(const int n, const int l) const;
};

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_STENCIL_ORTHOGONALVALUE_HPP
