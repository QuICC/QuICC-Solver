
#include <complex>
#include <iostream>

#include "Impl.hpp"
#include "View/View.hpp"
#include "ViewOps/Blas/Cpu/Gemm.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/TypeTraits.hpp"
#include "ViewOps/Worland/Types.hpp"

namespace QuICC {
namespace Transform {
namespace Worland {
namespace Cpu {

template <class Tout, class Tin, class Top, std::uint16_t Treatment>
void ImplOp<Tout, Tin, Top, Treatment>::applyImpl(Tout& out, const Tin& in,
   const Top& op)
{
   assert(op.dims()[0] == out.dims()[0]);
   assert(op.dims()[1] == in.dims()[0]);
   assert(op.dims()[2] == in.dims()[2]);
   assert(op.dims()[2] == out.dims()[2]);

   using IndexType = typename Tin::IndexType;


   ViewBase<IndexType> modsPointers;
   using LevelType = typename Top::LevelType;
   if constexpr (std::is_same_v<LevelType, CSL3D::level>)
   {
      // uniform truncation projector/integrator
      modsPointers = in.pointers()[1];
   }
   else
   {
      throw std::logic_error("backend for these types is not implemented.");
   }

   // cache populated layers
   if (_harmOrd.size() < 1)
   {
      for (IndexType k = 0; k < modsPointers.size() - 1; ++k)
      {
         IndexType nCols = modsPointers[k + 1] - modsPointers[k];
         // check if layer is populated
         if (nCols > 0)
         {
            _harmOrd.push_back(k);
            _cols.push_back(nCols);
         }
      }
   }

   // matmul loop
   std::uint32_t offSetA = 0;
   std::uint32_t offSetB = 0;
   std::uint32_t offSetC = 0;

   for (IndexType h = 0; h < _harmOrd.size(); ++h)
   {
      // select correct type to extract slice
      constexpr bool isSliceOpRowMaj =
         std::is_same_v<typename Top::OrderType, LoopOrderType<j_t, i_t, k_t>>;
      using opSliceAtt_t =
         std::conditional_t<isSliceOpRowMaj, dense2DRM, dense2D>;
      using dataSliceAtt_t =
         std::conditional_t<isSliceOpRowMaj, dense2D, dense2DRM>;

      IndexType M, N, K;

      // get dimensions
      if constexpr (std::is_same_v<LevelType, CSL3D::level>)
      {
         // uniform truncation projector/integrator
         M = op.dims()[0]; // radial points/modes
         K = op.dims()[1]; // radias modes/points
         N = _cols[h];     // number of columns
      }

      // check mem bounds
      assert(offSetA + M * K <= op.size());
      assert(offSetB + K * N <= in.size());
      assert(offSetC + M * N <= out.size());

      // set dense views
      View<typename Top::ScalarType, opSliceAtt_t> A(
         {op.data() + offSetA, M * K}, {M, K});
      View<typename Tin::ScalarType, dataSliceAtt_t> B(
         {in.data() + offSetB, K * N}, {K, N});
      View<typename Tout::ScalarType, dataSliceAtt_t> C(
         {out.data() + offSetC, M * N}, {M, N});

      // compute
      // if constexpr (Treatment == none_m)
      // {
      QuICC::Blas::Cpu::Eigen::matmul(C, A, B, 1.0);
      // }
      // else if constexpr (Treatment == diffPhi_m)
      // {
      //     std::complex<double> alpha{0.0, static_cast<double>(_harmOrd[h])};
      //     if constexpr (is_integrator_v<Top>)
      //     {
      //         alpha = -alpha;
      //     }
      //     QuICC::Blas::Cpu::Eigen::matmul(C, A, B, alpha);
      // }

      // update offset
      offSetA += A.size();
      offSetB += B.size();
      offSetC += C.size();
   }
}

// Explicit instantations

// Projectors and integrators (same data layout) with column major data
template class ImplOp<Uniform::phys_t, Uniform::mods_t, Uniform::projRM_t>;

// Projectors and integrators (same data layout) with row major data
template class ImplOp<Uniform::physRM_t, Uniform::modsRM_t, Uniform::proj_t>;


} // namespace Cpu
} // namespace Worland
} // namespace Transform
} // namespace QuICC
