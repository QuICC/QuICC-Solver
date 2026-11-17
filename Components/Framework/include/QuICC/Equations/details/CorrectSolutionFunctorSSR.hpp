/**
 * @file CorrectSolutionFunctorSSR.hpp
 * @brief Implementation of the CorrectSolution functor fpr SLOWEST_SINGLE_RHS
 */

#ifndef QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTORSSR_HPP
#define QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTORSSR_HPP

// System includes
//


// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/details/CorrectSolutionFunctor.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData>
void CorrectSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(
   TData& storage, const int start)
{
   for (auto&& c: corr)
   {
      Arithmetics::assignScalar<Arithmetics::Operation::Plus>(storage,
            std::get<1>(c) + start, std::get<0>(c));
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTORSSR_HPP
