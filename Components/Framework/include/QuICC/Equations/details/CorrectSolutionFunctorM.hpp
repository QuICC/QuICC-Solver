/**
 * @file CorrectSolutionFunctorM.hpp
 * @brief Implementation of the CorrectSolution functor fpr MODE
 */

#ifndef QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTORM_HPP
#define QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTORM_HPP

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
void CorrectSolutionFunctor<CouplingIndexType::MODE>::apply(TData& storage,
   const int start)
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

#endif // QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTORM_HPP
