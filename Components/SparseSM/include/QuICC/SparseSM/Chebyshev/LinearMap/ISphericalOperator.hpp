/**
 * @file ISphericalOperator.hpp
 * @brief Base implementation of a spherical operator, with y = ax + b and
 * needing the harmonic degree l
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_ISPHERICALOPERATOR_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_ISPHERICALOPERATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Base implementation of a spherical operator, with y = ax + b and
 * needing the harmonic degree l
 */
class ISphericalOperator : public ILinearMapOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param lower   Lower bound
    * @param upper   Upper bound
    * @param l       Harmonic degree l
    */
   ISphericalOperator(const int rows, const int cols, const Scalar_t lower,
      const Scalar_t upper, const Scalar_t l);

   /**
    * @brief Destructor
    */
   virtual ~ISphericalOperator() = default;

protected:
   /**
    * @brief Harmonic degree
    */
   Scalar_t l() const;

   /**
    * @brief Templated power function l^p
    */
   template <int P> Scalar_t l() const;

private:
   /**
    * @brief Harmonic degree l
    */
   Scalar_t mL;
};

template <int P> ISphericalOperator::Scalar_t ISphericalOperator::l() const
{
   if constexpr (P == 1)
   {
      return this->mL;
   }
   else if constexpr (P == 2)
   {
      return this->mL * this->mL;
   }
   else
   {
      return Internal::Math::pow(this->mL, P);
   }
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_ISPHERICALOPERATOR_HPP
