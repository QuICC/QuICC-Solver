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
#include "Arithmetics/Utility.hpp"

namespace QuICC {

namespace Equations {

class IFieldEquation;

namespace details
{
   template <typename TData>
   struct Temporary
   {
      TData data;
   };
}


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
   template <typename TData> Arithmetics::Temporary<TData> init(int& solStart, const TData& storage, const int start, const CouplingInformation& info);

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
template <typename TData> Arithmetics::Temporary<TData> StoreSolutionFunctor<IndexType>::init(int& solStart, const TData& storage, const int start, const CouplingInformation& info)
{
   Arithmetics::Temporary<TData> sol;
   if(info.isGalerkin())
   {
      if constexpr (Arithmetics::is_view<TData>::value)
      {
         sol.storage.resize(info.tauN(matIdx), info.rhsCols(matIdx));
         std::array<std::uint32_t, 2> dimensions {info.tauN(matIdx), info.rhsCols(matIdx)};
         Patch::std::span<typename TData::ScalarType> span(sol.storage.data(), sol.storage.size());
         sol.data = TData(span, storage.dims());
      }
      else
      {
         // Temporary storage is required
         sol.data = TData(info.tauN(matIdx), info.rhsCols(matIdx));
      }

      solStart = 0;
      sol.ptr = &sol.data;

      // Apply Galerkin stencil
      applyGalerkinStencil(*eq, compId, sol.data, start, matIdx, storage);
   }
   else
   {
      solStart = start;
      sol.ptr = &storage;
   }

   return sol;
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_STORESOLUTIONFUNCTOR_HPP
