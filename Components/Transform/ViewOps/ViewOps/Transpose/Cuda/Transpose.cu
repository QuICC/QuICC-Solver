/**
 * @file Transpose.cu
 * @brief Transpose operations on Views
 */

// External includes
//
#include <cassert>
#include <complex>

// Project includes
//
#include "Cuda/CudaUtil.hpp"
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "Transpose.hpp"
#include "View/View.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {

using namespace QuICC::Operator;



template <class Tout, class Tin, class Perm>
void Op<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("Reduction::Cuda::applyImpl");

    assert(QuICC::Cuda::isDeviceMemory(out.data()));
    assert(QuICC::Cuda::isDeviceMemory(in.data()));

}

// Explicit instantiations
template class Op<View::View<std::complex<double>, View::DCCSC3D>, View::View<std::complex<double>, View::DCCSC3DJIK>, p201_t>;


} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
