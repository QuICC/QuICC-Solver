/**
 * @file OrthoCgs2Functor.hpp
 * @brief Classical Gram-Schmidt orthogonalization with re-orthogonalization
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGS2FUNCTOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGS2FUNCTOR_HPP

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
 * @brief Classical Gram-Schmidt with re-orthogonalization
 */
template <typename TInner = StandardInnerProduct>
class OrthoCgs2Functor
{
public:
   /// @brief ctor
   OrthoCgs2Functor(std::shared_ptr<TInner> pInner, const int p);

   /// @brief ctor with Standard inner product
   template <typename T1 = TInner, typename  = typename std::enable_if_t<std::is_same_v<T1, StandardInnerProduct>>>
   OrthoCgs2Functor(const int p);

   /// @brief dtor
   ~OrthoCgs2Functor() = default;

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
OrthoCgs2Functor<TInner>::OrthoCgs2Functor(std::shared_ptr<TInner> pInner, const int p)
   : mcP(p), mpInner(pInner)
{}

template <typename TInner>
template <typename, typename>
OrthoCgs2Functor<TInner>::OrthoCgs2Functor(const int p)
   : mcP(p)
{
   this->mpInner = std::make_shared<StandardInnerProduct>();
}

template <typename TInner>
double OrthoCgs2Functor<TInner>::operator()(Matrix& matV, Matrix& matH, const int j, const int n)
{
   auto&& inner_product = *this->mpInner;

   // CGS Orthogonalization
   int i0 = std::max(0, j + 1 - this->mcP);
   Array colH = inner_product(matV, i0, j, matV, j+1, n);
   for(int i = i0; i <= j; i++)
   {
      matH(i, j) = colH(i-i0);

      matV.col(j+1) -= matH(i,j)*matV.col(i);
   }

   // CGS re-orthogonalization
   colH = inner_product(matV, i0, j, matV, j+1, n);
   for(int i = i0; i <= j; i++)
   {
      matH(i, j) += colH(i-i0);

      matV.col(j+1) -= colH(i-i0, 0)*matV.col(i);
   }

   // Norm
   auto normV = inner_product.norm(matV, j+1, n);
   return normV;
}

template <typename TInner>
int OrthoCgs2Functor<TInner>::p() const
{
   return this->mcP;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGS2FUNCTOR_HPP
