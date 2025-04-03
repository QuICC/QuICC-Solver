/**
 * @file Op.cu
 * @brief Transpose operations on Views
 */

// External includes
//
#include <cassert>
#include <complex>

// Project includes
//
#include "Cuda/CudaUtil.hpp"
#include "Op.hpp"
#include "Impl.hpp"
#include "Profiler/Interface.hpp"
#include "View/View.hpp"

#define QUICC_MAX_TH_NAIVE 2048

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {

using namespace QuICC::Operator;

template <class Tout, class Tin, class Perm>
void Op<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Transpose::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(out.data()));
   assert(QuICC::Cuda::isDeviceMemory(in.data()));

   if constexpr (std::is_same_v<Perm, p201_t>)
   {
      details::implPerm201(out, in);
   }
   else if constexpr (std::is_same_v<Perm, p120_t>)
   {
      details::implPerm120(out, in);
   }
   else
   {
      throw std::logic_error("transpose not implemented");
   }
}


// Explicit instantiations
// FT -> AL
template class Op<View::View<double, View::DCCSC3DJIK>,
   View::View<double, View::DCCSC3D>, p201_t>;
template class Op<View::View<std::complex<double>, View::DCCSC3DJIK>,
   View::View<std::complex<double>, View::DCCSC3D>, p201_t>;
// AL -> FT
template class Op<View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3DJIK>, p120_t>;
template class Op<View::View<std::complex<double>, View::DCCSC3D>,
   View::View<std::complex<double>, View::DCCSC3DJIK>, p120_t>;

// AL -> JW
template class Op<View::View<double, View::DCCSC3DJIK>,
   View::View<double, View::S1CLCSC3DJIK>, p201_t>;
template class Op<View::View<std::complex<double>, View::DCCSC3DJIK>,
   View::View<std::complex<double>, View::S1CLCSC3DJIK>, p201_t>;
// JW -> AL
template class Op<View::View<double, View::S1CLCSC3DJIK>,
   View::View<double, View::DCCSC3DJIK>, p120_t>;
template class Op<View::View<std::complex<double>, View::S1CLCSC3DJIK>,
   View::View<std::complex<double>, View::DCCSC3DJIK>, p120_t>;

} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
