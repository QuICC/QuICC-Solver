/**
 * @file StoreSolutionFunctor.hpp
 * @brief Implementation of the CopyUnknown functors
 */

#ifndef QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/Utility.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/ApplyGalerkinStencil.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/SolutionUpdater.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class StoreSolutionFunctor
{
public:
   /**
    * @brief ctor
    */
   StoreSolutionFunctor(const Resolution& res, const CouplingInformation& cinfo, const SparseMatrix* op, std::shared_ptr<SolutionUpdater> spUp,
      const int matIdx);

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
    * @param storage Solver solution
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField>
   void apply(TField& field, const TData& storage, const int start);

private:
   template <typename TData>
   Arithmetics::Temporary<TData> init(int& solStart, const TData& storage,
      const int start);

   /**
    * @brief Resolution
    */
   const Resolution& res;

   /**
    * @brief Coupling information
    */
   const CouplingInformation& cinfo;

   /**
    * @brief Reference to equation
    */
   const SparseMatrix* pOp;

   /**
    * @brief Update functor
    */
   std::shared_ptr<SolutionUpdater> spUp;

   /**
    * @brief Matrix index
    */
   const int matIdx;
};

template <CouplingIndexType IndexType>
StoreSolutionFunctor<IndexType>::StoreSolutionFunctor(const Resolution& res, const CouplingInformation& cinfo, const SparseMatrix* pOp, std::shared_ptr<SolutionUpdater> spUp,
   const int matIdx) :
    res(res), cinfo(cinfo), pOp(pOp), spUp(spUp), matIdx(matIdx)
{}

template <CouplingIndexType IndexType>
template <typename TData>
Arithmetics::Temporary<TData> StoreSolutionFunctor<IndexType>::init(
   int& solStart, const TData& storage, const int start)
{
   Arithmetics::Temporary<TData> sol;
   if (cinfo.isGalerkin())
   {
      if constexpr (Arithmetics::is_view<TData>::value)
      {
         sol.storage.resize(cinfo.tauN(matIdx), cinfo.rhsCols(matIdx));
         std::array<std::uint32_t, 2> dimensions{cinfo.tauN(matIdx),
            cinfo.rhsCols(matIdx)};
         Patch::std::span<typename TData::ScalarType> span(sol.storage.data(),
            sol.storage.size());
         sol.data = TData(span, storage.dims());
      }
      else
      {
         // Temporary storage is required
         sol.data = TData(cinfo.tauN(matIdx), cinfo.rhsCols(matIdx));
      }

      solStart = 0;
      sol.ptr = &sol.data;

      // Apply Galerkin stencil
      applyGalerkinStencil(*pOp, sol.data, start, matIdx, storage);
   }
   else
   {
      solStart = start;
      sol.ptr = &storage;
   }

   return sol;
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTOR_HPP
