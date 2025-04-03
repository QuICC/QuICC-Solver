/**
 * @file OpGrouped.cu
 * @brief Grouped Transpose operations on Views
 */

// External includes
//
#include <cassert>
#include <complex>

// Project includes
//
#include "Cuda/CudaUtil.hpp"
#include "OpGrouped.hpp"
#include "Impl.hpp"
#include "Profiler/Interface.hpp"
#include "View/View.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {


template <class Tout, class Tin, class Perm>
void OpGrouped<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
   assert(out.size() == in.size());
   Profiler::RegionFixture<4> fix("Transpose::Cuda::OpGrouped::applyImpl");

   if constexpr (std::is_same_v<Perm, p201_t>)
   {
      for (std::size_t i = 0; i < in.size(); ++i)
      {
         assert(QuICC::Cuda::isDeviceMemory(out[i].data()));
         assert(QuICC::Cuda::isDeviceMemory(in[i].data()));
         details::implPerm201(out[i], in[i]);
      }
   }
   else if constexpr (std::is_same_v<Perm, p120_t>)
   {
      for (std::size_t i = 0; i < in.size(); ++i)
      {
         assert(QuICC::Cuda::isDeviceMemory(out[i].data()));
         assert(QuICC::Cuda::isDeviceMemory(in[i].data()));
         details::implPerm120(out[i], in[i]);
      }
   }
   else
   {
      throw std::logic_error("transpose not implemented");
   }
}


// Explicit instantiations
// FT -> AL
template class OpGrouped<std::vector<View::View<double, View::DCCSC3DJIK>>,
   std::vector<View::View<double, View::DCCSC3D>>, p201_t>;
template class OpGrouped<std::vector<View::View<std::complex<double>, View::DCCSC3DJIK>>,
   std::vector<View::View<std::complex<double>, View::DCCSC3D>>, p201_t>;
// AL -> FT
template class OpGrouped<std::vector<View::View<double, View::DCCSC3D>>,
   std::vector<View::View<double, View::DCCSC3DJIK>>, p120_t>;
template class OpGrouped<std::vector<View::View<std::complex<double>, View::DCCSC3D>>,
   std::vector<View::View<std::complex<double>, View::DCCSC3DJIK>>, p120_t>;

// AL -> JW
template class OpGrouped<std::vector<View::View<double, View::DCCSC3DJIK>>,
   std::vector<View::View<double, View::S1CLCSC3DJIK>>, p201_t>;
template class OpGrouped<std::vector<View::View<std::complex<double>, View::DCCSC3DJIK>>,
   std::vector<View::View<std::complex<double>, View::S1CLCSC3DJIK>>, p201_t>;
// JW -> AL
template class OpGrouped<std::vector<View::View<double, View::S1CLCSC3DJIK>>,
   std::vector<View::View<double, View::DCCSC3DJIK>>, p120_t>;
template class OpGrouped<std::vector<View::View<std::complex<double>, View::S1CLCSC3DJIK>>,
   std::vector<View::View<std::complex<double>, View::DCCSC3DJIK>>, p120_t>;

} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
