#include <iostream>
#include <complex>

#include "Spec.hpp"
#include "View/View.hpp"
#include "ViewOps/Chebyshev/LinearMap/Util.hpp"
#include "ViewOps/Chebyshev/LinearMap/Tags.hpp"
#include "ViewOps/Chebyshev/LinearMap/Types.hpp"
#include "Profiler/Interface.hpp"

#ifdef QUICC_HAS_CUDA_BACKEND
#include "Cuda/CudaUtil.hpp"
#endif

namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {
namespace Cpu {

template<class Tout, class Tin, class Operation, std::uint16_t Treatment>
SpecOp<Tout, Tin, Operation, Treatment>::SpecOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, class Operation, std::uint16_t Treatment>
void SpecOp<Tout, Tin, Operation, Treatment>::applyImpl(Tout& out, const Tin& in, const ScaleType fftScaling)
{
    Profiler::RegionFixture<5> fix("SpecOp::applyImpl");

    assert(out.size() == in.size());
    assert(out.dims()[0] == in.dims()[0]);
    assert(out.dims()[1] == in.dims()[1]);
    assert(out.dims()[2] == in.dims()[2]);

    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(!QuICC::Cuda::isDeviceMemory(out.data()));
    #endif

    if constexpr (std::is_same_v<Operation, spec_id>)
    {
        // if the spec is identity and in place and there are no modes
        // to be zeroed then it is a noop
        if(out.data() == in.data() && out.dims()[0] == out.lds())
        {
            return;
        }
    }

    float c = fftScaling;

    // Column major
    // Get total number of columns to loop over
    auto indices = in.indices()[1];
    auto columns = indices.size();

    const auto N = in.lds();
    const auto nDealias = in.dims()[0];
    for (std::size_t col = 0; col < columns ; ++col)
    {
        // linear index (:,n,k)
        std::size_t nk = N*col;
        std::size_t n = 0;

        for (; n < nDealias; ++n)
        {
           out.data()[nk+n] = in.data()[nk+n];
        }

        for (; n < N; ++n)
        {
            out.data()[nk+n] = 0.0;
        }

    }

}

// explicit instantations
template class SpecOp<mods_t, mods_t, spec_id>;

} // namespace Cpu
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
