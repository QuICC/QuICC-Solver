#include <complex>
#include <iostream>

#include "Grid.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"
#include "View/View.hpp"
#include "ViewOps/Chebyshev/LinearMap/Tags.hpp"
#include "ViewOps/Chebyshev/LinearMap/Types.hpp"
#include "ViewOps/Chebyshev/LinearMap/Util.hpp"

#ifdef QUICC_HAS_CUDA_BACKEND
#include "Cuda/CudaUtil.hpp"
#endif

namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {
namespace Cpu {

template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
GridOp<Tout, Tin, Operation, Treatment>::GridOp(const double lower,
   const double upper) :
    mLower(lower), mUpper(upper), mScale(1.0){};

template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
GridOp<Tout, Tin, Operation, Treatment>::GridOp(const double lower,
   const double upper, ScaleType scale) :
    mLower(lower), mUpper(upper), mScale(scale){};

template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
void GridOp<Tout, Tin, Operation, Treatment>::applyImpl(Tout& out,
   const Tin& in, const ScaleType fftScaling)
{
   Profiler::RegionFixture<5> fix("GridOp::applyImpl");

   assert(out.size() == in.size());
   assert(out.dims()[0] == in.dims()[0]);
   assert(out.dims()[1] == in.dims()[1]);
   assert(out.dims()[2] == in.dims()[2]);

#ifdef QUICC_HAS_CUDA_BACKEND
   assert(!QuICC::Cuda::isDeviceMemory(out.data()));
#endif

   if constexpr (std::is_same_v<Operation, grid_id>)
   {
      // if the grid is identity and in place and there are no modes
      // to be zeroed then it is a noop
      if (out.data() == in.data() && out.dims()[0] == out.lds())
      {
         return;
      }
   }

   // Column major
   // Get total number of columns to loop over
   auto indices = in.indices()[1];
   auto columns = indices.size();

   const auto N = in.lds();
   const auto nDealias = in.dims()[0];
   if constexpr (std::is_same_v<Operation, grid_divy1> ||
                 std::is_same_v<Operation, grid_divy2>)
   {
      // Initialize grid scaling
      if (mGridScaler.size() == 0)
      {
         int p = 0;
         if constexpr (std::is_same_v<Operation, grid_divy1>)
         {
            p = -1;
         }
         else if constexpr (std::is_same_v<Operation, grid_divy2>)
         {
            p = -2;
         }

         Internal::Array igrid, iweights;
         Polynomial::Quadrature::ChebyshevRule quad;
         quad.computeQuadrature(igrid, iweights, nDealias, mLower, mUpper);

         mGridScaler.reserve(nDealias);
         for (std::size_t i = 0; i < nDealias; i++)
         {
            mGridScaler.push_back(
               static_cast<double>(Internal::Math::pow(igrid(i), p)));
         }
      }

      for (std::size_t col = 0; col < columns; ++col)
      {
         // linear index (:,n,k)
         std::size_t nk = N * col;
         std::size_t n = 0;

         for (; n < nDealias; ++n)
         {
            out.data()[nk + n] = in.data()[nk + n] * mGridScaler[n];
         }
      }
   }
   else
   {
      for (std::size_t col = 0; col < columns; ++col)
      {
         // linear index (:,n,k)
         std::size_t nk = N * col;
         std::size_t n = 0;

         for (; n < nDealias; ++n)
         {
            out.data()[nk + n] = in.data()[nk + n];
         }
      }
   }
}

// explicit instantations
template class GridOp<phys_t, phys_t, grid_id>;
template class GridOp<phys_t, phys_t, grid_divy1>;
template class GridOp<phys_t, phys_t, grid_divy2>;

} // namespace Cpu
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
