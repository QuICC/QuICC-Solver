/**
 * @file ILinearMapOperator.hpp
 * @brief Implementation of the generic interface to the Cheyshev sparse
 * operator based on a linear map y = ax + b, x = [-1, 1] (natural chebyshev
 * grid)
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_ILINEARMAPOPERATOR_HPP
#define QUICC_SPARSESM_CHEBYSHEV_ILINEARMAPOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/SparseSM/IChebyshevOperator.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

/**
 * @brief Implementation of the generic interface to the Chebyshev sparse
 * operator based on a linear map y = ax + b, x = [-1, 1] (natural chebyshev
 * grid)
 */
class ILinearMapOperator : public IChebyshevOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of columns
    * @param lower   Lower bound of y
    * @param upper   Lower bound of y
    */
   ILinearMapOperator(const int rows, const int cols, const Scalar_t lower,
      const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~ILinearMapOperator();

protected:
   /**
    * @brief Get mapping a coefficient from y = ax + b
    */
   Scalar_t a() const;

   /**
    * @brief Templated power function a^p
    */
   template <int P> Scalar_t a() const;

   /**
    * @brief Get mapping b coefficient from y = ax + b
    */
   Scalar_t b() const;

   /**
    * @brief Templated power function b^p
    */
   template <int P> Scalar_t b() const;

private:
   /**
    * @brief Compute mapping
    *
    * @param lower   Lower bound
    * @param upper   Upper bound
    */
   void setBounds(const Scalar_t lower, const Scalar_t upper);

   /**
    * @brief a coefficienct of y = ax + b
    */
   Scalar_t mA;

   /**
    * @brief b coefficienct of y = ax + b
    */
   Scalar_t mB;
};

template <int P> ILinearMapOperator::Scalar_t ILinearMapOperator::a() const
{
   if constexpr (P == 1)
   {
      return this->mA;
   }
   else if constexpr (P == 2)
   {
      return this->mA * this->mA;
   }
   else
   {
      return Internal::Math::pow(this->mA, P);
   }
}

template <int P> ILinearMapOperator::Scalar_t ILinearMapOperator::b() const
{
   if constexpr (P == 1)
   {
      return this->mB;
   }
   else if constexpr (P == 2)
   {
      return this->mB * this->mB;
   }
   else
   {
      return Internal::Math::pow(this->mB, P);
   }
}

} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_ILINEARMAPOPERATOR_HPP
