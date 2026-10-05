/**
 * @file OrthoMgsFunctor.hpp
 * @brief Modified Gram-Schmidt orthogonalization
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOMGSFUNCTOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOMGSFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Timestep/Exponential/Functors/StandardInnerProduct.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

namespace Functors {

/**
 * @brief Modified Gram-Schmidt
 */
template <typename TInner = StandardInnerProduct>
class OrthoMgsFunctor
{
public:
   /// @brief ctor
   OrthoMgsFunctor(std::shared_ptr<TInner> pInner, const int p);

   /// @brief ctor with Standard inner product
   template <typename T1 = TInner, typename  = typename std::enable_if_t<std::is_same_v<T1, StandardInnerProduct>>>
   OrthoMgsFunctor(const int p);

   /// @brief dtor
   ~OrthoMgsFunctor() = default;

   /**
    * Orthogonalize
    */
   double operator()(Matrix& matV, Matrix& matH, const int j, const int n);

   /**
    * @brief Gram-schmidt order
    */
   int p() const;

private:
   /**
    * @brief Length of incomplete orthogonalization
    */
   const int mcP;

   /**
    * @brief Inner product functor
    */
   std::shared_ptr<TInner> mpInner;
};

template <typename TInner>
OrthoMgsFunctor<TInner>::OrthoMgsFunctor(std::shared_ptr<TInner> pInner, const int p)
   : mcP(p), mpInner(pInner)
{}

template <typename TInner>
template <typename, typename>
OrthoMgsFunctor<TInner>::OrthoMgsFunctor(const int p)
   : mcP(p)
{
   this->mpInner = std::make_shared<StandardInnerProduct>();
}

template <typename TInner>
double OrthoMgsFunctor<TInner>::operator()(Matrix& matV, Matrix& matH, const int j, const int n)
{
   auto&& inner_product = *this->mpInner;

   // MGS Orthogonalization
   int i0 = std::max(0, j + 1 - this->mcP);
   for(int i = i0; i <= j; i++)
   {
      auto hij = inner_product(matV, i, matV, j+1, n);

      matH(i, j) = hij;

      matV.col(j+1) -= hij*matV.col(i);
   }

   // Norm
   auto normV = inner_product.norm(matV, j+1, n);
   return normV;
}

template <typename TInner>
int OrthoMgsFunctor<TInner>::p() const
{
   return this->mcP;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOMGSFUNCTOR_HPP
