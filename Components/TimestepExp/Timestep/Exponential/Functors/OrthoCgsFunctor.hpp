/**
 * @file OrthoCgsFunctor.hpp
 * @brief Classical Gram-Schmidt orthogonalization
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGSFUNCTOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGSFUNCTOR_HPP

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
 * @brief Classical Gram-Schmidt
 */
template <typename TInner = StandardInnerProduct>
class OrthoCgsFunctor
{
public:
   /// @brief ctor
   OrthoCgsFunctor(std::shared_ptr<TInner> pInner, const int p);

   /// @brief ctor with Standard inner product
   template <typename T1 = TInner, typename  = typename std::enable_if_t<std::is_same_v<T1, StandardInnerProduct>>>
   OrthoCgsFunctor(const int p);

   /// @brief dtor
   ~OrthoCgsFunctor() = default;

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
OrthoCgsFunctor<TInner>::OrthoCgsFunctor(std::shared_ptr<TInner> pInner, const int p)
   : mcP(p), mpInner(pInner)
{}

template <typename TInner>
template <typename, typename>
OrthoCgsFunctor<TInner>::OrthoCgsFunctor(const int p)
   : mcP(p)
{
   this->mpInner = std::make_shared<StandardInnerProduct>();
}

template <typename TInner>
double OrthoCgsFunctor<TInner>::operator()(Matrix& matV, Matrix& matH, const int j, const int n)
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

   // Norm
   auto normV = inner_product.norm(matV, j+1, n);
   return normV;
}

template <typename TInner>
int OrthoCgsFunctor<TInner>::p() const
{
   return this->mcP;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGSFUNCTOR_HPP
