/**
 * @file StoreSolutionFunctor.hpp
 * @brief Implementation of the CopyUnknown functors
 */

#ifndef QUICC_EQUATIONS_STORESOLUTIONFUNCTOR_HPP
#define QUICC_EQUATIONS_STORESOLUTIONFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"

namespace QuICC {

namespace Equations {

class IFieldEquation;

template <CouplingIndexType IndexType>
class StoreSolutionFunctor
{
 public:
   /**
    * @brief ctor
    */
   StoreSolutionFunctor(IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief deleted default ctor
    */
   StoreSolutionFunctor() = delete;

   /**
    * @brief dtor
    */
   ~StoreSolutionFunctor() = default;

   /**
    * @brief Transfer solver solution to equation unknown
    *
    * @param field   Scalar or vector field
    * @param compId  Component ID
    * @param storage Solver solution
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void apply(TField& field, const TData& storage, const int start);

 private:

   template <typename TData> const TData* init(int& solStart, const TData& storage, const int start, TData& tmp, const CouplingInformation& info);

   /**
    * @brief Reference to equation
    */
   IFieldEquation* eq;

   /**
    * @brief Field component ID
    */
   FieldComponents::Spectral::Id compId;

   /**
    * @brief Matrix index
    */
   const int matIdx;
};

template <CouplingIndexType IndexType>
StoreSolutionFunctor<IndexType>::StoreSolutionFunctor(IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   : eq(&eq), compId(compId), matIdx(matIdx)
{
}

template <CouplingIndexType IndexType>
template <typename TData> const TData* StoreSolutionFunctor<IndexType>::init(int& solStart, const TData& storage, const int start, TData& tmp, const CouplingInformation& info)
{
   const TData * solution;
   if(info.isGalerkin())
   {
      // Temporary storage is required
      tmp = TData(info.tauN(matIdx), info.rhsCols(matIdx));

      // Apply Galerkin stencil
      applyGalerkinStencil(*eq, compId, tmp, start, matIdx, storage);

      solStart = 0;
      solution = &tmp;
   }
   else
   {
      solStart = start;
      solution = &storage;
   }

   return solution;
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_STORESOLUTIONFUNCTOR_HPP
