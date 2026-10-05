/**
 * @file StandardInnerProduct.hpp
 * @brief Functor for energy based innner product
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_STANDARDINNERPRODUCT_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_STANDARDINNERPRODUCT_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Functor for energy based inner product
 */
class StandardInnerProduct
{
   public:
      /**
       * @brief ctor
       */
      StandardInnerProduct() = default;

      /**
       * @brief ctor
       */
      ~StandardInnerProduct() = default;

      /**
       * @brief Compute inner product
       */
      MHDFloat operator()(const Matrix& u, const int i, const Matrix& v, const int j, const int n)  const;

      /**
       * @brief Compute inner product for multiple vectors
       */
      Array operator()(const Matrix& u, const int i0, const int i1, const Matrix& v, const int j, const int n)  const;

      /**
       * @brief Compute norm based on innner product
       */
      MHDFloat norm(const Matrix& u, const int i, const int n)  const;

      /**
       * @brief Set data handle
       */
      void setWorkspaceHandle(Matrix& tmp) {};

   private:
};

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_STANDARDINNERPRODUCT_HPP
