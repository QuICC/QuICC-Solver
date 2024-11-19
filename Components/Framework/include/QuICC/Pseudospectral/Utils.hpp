/**
 * @file Utils.hpp
 * @brief Pseudospectral Utils
 */

#ifndef QUICC_PSEUDOSPECTRAL_UTILS_HPP
#define QUICC_PSEUDOSPECTRAL_UTILS_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "View/View.hpp"
#include "Memory/Memory.hpp"
#include "QuICC/Resolutions/TransformResolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "Graph/Types.hpp"

namespace QuICC {

namespace Pseudospectral {

std::size_t hash_combine(const std::size_t a, const std::size_t b);

namespace details
{
   struct ptrAndIdxBlock
   {
      Memory::MemBlock<std::uint32_t> ptr;
      Memory::MemBlock<std::uint32_t> idx;
   };

   ptrAndIdxBlock getMeta(const TransformResolution& res, const std::uint32_t maxLayers, std::shared_ptr<Memory::memory_resource> mem);

   void copyScalar2View(Graph::varData_t vVar, const Framework::Selector::VariantSharedScalarVariable sVar, const TransformResolution& res);

   void copyVector2View(Graph::varData_t vVar0, Graph::varData_t vVar1, const Framework::Selector::VariantSharedVectorVariable sVar, const TransformResolution& res);

   void copyView2Scalar(Framework::Selector::VariantSharedScalarVariable sVar, const Graph::varData_t vecVar, const TransformResolution& res);

   void copyView2Vector(Framework::Selector::VariantSharedVectorVariable vecVar, const Graph::varData_t vVar0, Graph::varData_t vVar1, const TransformResolution& res);

   void copyView2Vector(Framework::Selector::VariantSharedVectorVariable vecVar, const Graph::varData_t vVar0, const Graph::varData_t vVar1, const Graph::varData_t vVar2, const TransformResolution& res);

} // namespace details
} // Pseudospectral
} // QuICC

#endif // QUICC_PSEUDOSPECTRAL_UTILS_HPP
