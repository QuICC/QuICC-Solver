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
    * @param dBeta   Jacobi dBeta parameter
    * @param l       Harmonic degree l
    * @param nId     Normalization ID
    * @param c       Scaling constant
    */
   OrthogonalValue(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l, const std::size_t nId, const Scalar_t c = 1);

   /**
    * @brief Destructor
    */
   virtual ~OrthogonalValue() = default;

protected:
   enum class OrthoId: std::size_t {
      TorSphEnergy = 0,
      PolSphEnergy,
      ScaSphEnergy
   };

   /**
    * @brief Implementation of build dense matrix operator
    * @param output operator
    */
   virtual void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const override;

   /**
    * @brief Basis base Jacobi alpha
    */
   static Scalar_t basisAlpha(const std::size_t nId);

   /**
    * @brief Basis base Jacobi alpha
    */
   static Scalar_t basisDBeta(const std::size_t nId);

private:
   /**
    * @brief Normalization factor
    */
   Scalar_t norm(const int n, const int l) const;
};

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_STENCIL_ORTHOGONALVALUE_HPP
