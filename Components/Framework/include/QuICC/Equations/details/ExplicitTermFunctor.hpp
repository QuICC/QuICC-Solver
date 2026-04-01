/**
 * @file ExplicitTermFunctor.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class ExplicitTermFunctor
{
public:
   /**
    * @brief ctor
    */
   ExplicitTermFunctor(const Resolution& res, const CouplingInformation& cinfo,
      const int matIdx);

   /**
    * @brief deleted default ctor
    */
   ExplicitTermFunctor() = delete;

   /**
    * @brief dtor
    */
   ~ExplicitTermFunctor() = default;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param explicitField Explicit linear field values
    */
   template <typename T, typename TOperator, typename TData>
   void apply(TData& rSolverField, const TOperator& op, const int eqStart,
      const Framework::Selector::ScalarField<T>& explicitField);

private:
   /**
    * @brief Resolution
    */
   const Resolution& res;

   /**
    * @brief Coupling information
    */
   const CouplingInformation& cinfo;

   /**
    * @brief Matrix index
    */
   const int matIdx;
};

template <CouplingIndexType IndexType>
ExplicitTermFunctor<IndexType>::ExplicitTermFunctor(const Resolution& res, const CouplingInformation& cinfo,
   const int matIdx) :
    res(res), cinfo(cinfo), matIdx(matIdx)
{}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTOR_HPP
