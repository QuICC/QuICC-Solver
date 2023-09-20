
#include <complex>
#include <iostream>
#include <Eigen/Core>

#include "Impl.hpp"
#include "View/View.hpp"
#include "ViewOps/ALegendre/Tags.hpp"
#include "ViewOps/ALegendre/Types.hpp"
#include "ViewOps/ALegendre/TypeTraits.hpp"
#include "ViewOps/Blas/Cpu/Gemm.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {
namespace Cpu {

template<class Tout, class Tin, class Top, std::uint16_t Treatment>
void ImplOp<Tout, Tin, Top, Treatment>::applyImpl(Tout& out, const Tin& in, const Top& op)
{
    using IndexType = typename Tin::IndexType;

    IndexType L; // harmonic order

    ViewBase<IndexType> modsPointers;
    if constexpr (is_projector_v<Top>)
    {
        modsPointers = in.pointers()[1];
        L = in.dims()[0];
    }
    else if constexpr (is_integrator_v<Top>)
    {
        modsPointers = out.pointers()[1];
        L = out.dims()[0];
    }
    else
    {
        static_assert("backend for these types is not implemented.");
    }

    // cache populated layers
    if (_harmOrd.size() < 1)
    {
        for (IndexType k = 0; k < modsPointers.size()-1 ; ++k)
        {
            IndexType nCols = modsPointers[k+1] - modsPointers[k];
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

    for (IndexType h = 0; h < _harmOrd.size() ; ++h)
    {
        // select correct type to extract slice
        constexpr bool isSliceOpRowMaj = std::is_same_v<typename Top::OrderType, LoopOrderType<j_t, i_t, k_t>>;
        using opSliceAtt_t = std::conditional_t<isSliceOpRowMaj, dense2DRM, dense2D>;
        using dataSliceAtt_t = std::conditional_t<isSliceOpRowMaj, dense2D, dense2DRM>;

        IndexType M, N, K;

        // get dimensions
        if constexpr (is_projector_v<Top>)
        {
            M = out.dims()[0]; // longitudinal points
            K = L - _harmOrd[h]; // current harmonic order
            N = _cols[h]; // number of columns
        }
        else if constexpr (is_integrator_v<Top>)
        {
            M = L - _harmOrd[h]; // current harmonic order
            K = in.dims()[0]; // longitudinal points
            N = _cols[h]; // number of columns
        }
        else
        {
            static_assert("backend for these types is not implemented.");
        }

        // set dense views
        View<typename Top::ScalarType, opSliceAtt_t> A ({op.data() + offSetA, M*K}, {M, K});
        View<typename Tin::ScalarType, dataSliceAtt_t> B ({in.data() + offSetB, K*N},{K, N});
        View<typename Tout::ScalarType, dataSliceAtt_t> C ({out.data() + offSetC, M*N},{M, N});

        // compute
        if constexpr (Treatment == none_m)
        {
            QuICC::Blas::Cpu::Eigen::matmul(C, A, B, 1.0);
        }
        else if constexpr (Treatment == diffPhi_m)
        {
            std::complex<double> alpha{0.0, static_cast<double>(_harmOrd[h])};
            if constexpr (is_integrator_v<Top>)
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

// Explicit instantations
using namespace QuICC::Memory;

// Projectors with column major data
template class ImplOp<phys_t, mods_t, projRM_t>;
template class ImplOp<phys_t, mods_t, projRM_t, diffPhi_m>;

// Projectors with row major data
template class ImplOp<physRM_t, modsRM_t, proj_t>;
template class ImplOp<physRM_t, modsRM_t, proj_t, diffPhi_m>;

// Integrators with column major data
template class ImplOp<mods_t, phys_t, intRM_t>;
template class ImplOp<mods_t, phys_t, intRM_t, diffPhi_m>;

// Integrators with row major data
template class ImplOp<modsRM_t, physRM_t, int_t>;
template class ImplOp<modsRM_t, physRM_t, int_t, diffPhi_m>;

} // namespace Cpu
} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
