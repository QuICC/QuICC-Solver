/**
 * @file Impl.hpp
 * @brief cpu implementation of Associated Legendre operator
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Binary.hpp"
#include "View/View.hpp"
#include "ViewOps/Blas/Cpu/Gemm.hpp"
#include "ViewOps/Quadrature/Tags.hpp"
#include "ViewOps/Quadrature/ViewBatchedMatmulUtils.hpp"

namespace QuICC {
namespace Transform {
namespace Quadrature {
/// @brief Cpu backend namespace
namespace Cpu {

/// @brief Derived classes implement Associated Legendre operators.
/// In practice it boils down to a batched matmul with each batch having
/// different sizes.
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Top operator type
/// @tparam Treatment tag to include scaling due to derivative
template <class Tout, class Tin, class Top, std::uint16_t Treatment = 0>
class ImplOp
    : public Operator::BinaryBaseOp<ImplOp<Tout, Tin, Top, Treatment>, Tout, Tin, Top>
{
public:
   /// @brief Default constructor
   ImplOp() = default;
   /// @brief dtor
   ~ImplOp() = default;

private:
   /// @brief Action implementation
   /// @param out differentiatied modes
   /// @param in input modes
   /// @param op operator
   void applyImpl(Tout& out, const Tin& in, const Top& op);
   /// @brief Give access to base class
   friend Operator::BinaryBaseOp<ImplOp<Tout, Tin, Top, Treatment>, Tout, Tin, Top>;
   /// @brief index typedef
   using IndexType = typename Tin::IndexType;
   /// @brief layer index cache
   std::vector<IndexType> _layerIndex;
   /// @brief layer width cache
   std::vector<IndexType> _layerWidth;
};

template <class Tout, class Tin, class Top, std::uint16_t Treatment>
void ImplOp<Tout, Tin, Top, Treatment>::applyImpl(Tout& out, const Tin& in,
   const Top& op)
{
   // batched matmul out = op*in
   assert(op.dims()[0] == out.dims()[0]);
   assert(op.dims()[1] == in.dims()[0]);
   assert(op.dims()[2] == in.dims()[2]);
   assert(op.dims()[2] == out.dims()[2]);

   auto modsPointers = getModsPointers(out, in, op);

   using IndexType = typename Tin::IndexType;
   // cache populated layers
   if (_layerIndex.size() < 1)
   {
      for (IndexType k = 0; k < modsPointers.size() - 1; ++k)
      {
         IndexType nCols = modsPointers[k + 1] - modsPointers[k];
         assert(nCols <= in.dims()[1]);
         // check if layer is populated
         if (nCols > 0)
         {
            _layerIndex.push_back(k);
            _layerWidth.push_back(nCols);
         }
      }
   }

   // matmul loop
   std::uint32_t offSetA = 0;
   std::uint32_t offSetB = 0;
   std::uint32_t offSetC = 0;

   for (IndexType h = 0; h < _layerIndex.size(); ++h)
   {
      // select correct type to extract slice
      constexpr bool isSliceOpRowMaj =
         std::is_same_v<typename Top::OrderType, View::LoopOrderType<View::j_t, View::i_t, View::k_t>>;
      using opSliceAtt_t =
         std::conditional_t<isSliceOpRowMaj, View::dense2DRM, View::dense2D>;
      using dataSliceAtt_t =
         std::conditional_t<isSliceOpRowMaj, View::dense2D, View::dense2DRM>;

      // get dimensions
      auto dims = getMatmulDims(out, in, op, _layerWidth[h], _layerIndex[h]);
      auto M = dims.M;
      auto K = dims.K;
      auto N = dims.N;

      // check mem bounds
      assert(offSetA + M * K <= op.size());
      assert(offSetB + K * N <= in.size());
      assert(offSetC + M * N <= out.size());

      // set dense views
      View::View<typename Top::ScalarType, opSliceAtt_t> A(
         {op.data() + offSetA, M * K}, {M, K});
      View::View<typename Tin::ScalarType, dataSliceAtt_t> B(
         {in.data() + offSetB, K * N}, {K, N});
      View::View<typename Tout::ScalarType, dataSliceAtt_t> C(
         {out.data() + offSetC, M * N}, {M, N});

      // compute
      if constexpr (Treatment == none_m)
      {
         QuICC::Blas::Cpu::Eigen::matmul(C, A, B, 1.0);
      }
      else
      {
         std::complex<double> alpha{0.0, static_cast<double>(_layerIndex[h])};
         if constexpr (Treatment == diffPhiInt_m)
         {
            alpha = -alpha;
         }
         QuICC::Blas::Cpu::Eigen::matmul(C, A, B, alpha);
      }

      // update offset
      offSetA += A.size();
      offSetB += B.size();
      offSetC += C.size();
   }
}

} // namespace Cpu
} // namespace Quadrature
} // namespace Transform
} // namespace QuICC
