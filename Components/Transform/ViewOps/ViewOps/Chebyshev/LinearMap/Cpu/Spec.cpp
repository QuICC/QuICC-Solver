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
SpecOp<Tout, Tin, Operation, Treatment>::SpecOp(const double lower, const double upper) : mLower(lower), mUpper(upper), mScale(1.0){};

template<class Tout, class Tin, class Operation, std::uint16_t Treatment>
SpecOp<Tout, Tin, Operation, Treatment>::SpecOp(const double lower, const double upper, ScaleType scale) : mLower(lower), mUpper(upper), mScale(scale){};

template<class Tout, class Tin, class Operation, std::uint16_t Treatment>
void SpecOp<Tout, Tin, Operation, Treatment>::applyImpl(Tout& out, const Tin& in, const ScaleType fftScaling)
{
    Profiler::RegionFixture<5> fix("SpecOp::applyImpl");

    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(!QuICC::Cuda::isDeviceMemory(out.data()));
    #endif

    assert(out.dims()[1] == in.dims()[1]);
    assert(out.dims()[2] == in.dims()[2]);

    std::size_t nDealias;
    if constexpr (std::is_same_v<Operation, spec_id>)
    {
       assert(out.size() <= in.size());
       assert(out.dims()[0] <= in.dims()[0]);

       // if the spec is identity and in place and there are no modes
       // to be zeroed then it is a noop
       if(out.data() == in.data() && out.dims()[0] == out.lds())
       {
          return;
       }

       nDealias = out.dims()[0];
    }
    else if constexpr(std::is_same_v<Operation, spec_pad>)
    {
       assert(out.size() >= in.size());
       assert(out.dims()[0] >= in.dims()[0]);

       nDealias = in.dims()[0];
    }

    std::size_t Nout = out.lds();
    std::size_t Nin = in.lds();

    ScaleType c = fftScaling;

    // Column major
    // Get total number of columns to loop over
    auto indices = in.indices()[1];
    auto columns = indices.size();

    if constexpr(Operation::p.size() == 1 && Operation::p[0] == 0)
    {
       for (std::size_t col = 0; col < columns ; ++col)
       {
           // linear index (:,n,k)
           std::size_t nko = Nout*col;
           std::size_t nki = Nin*col;
           std::size_t n = 0;

           for (; n < nDealias; ++n)
           {
              out.data()[nko+n] = in.data()[nki+n] * c;
           }

           for (; n < Nout; ++n)
           {
              out.data()[nko+n] = 0;
           }
       }
    }
    else
    {
       const auto& ps = Operation::p;
       const auto& ts = Operation::t;

       Nout = out.lds();
       Nin = in.lds();
       nDealias = in.dims()[0];

       // 2*a, from y=ax + b
       double a2 = (mUpper - mLower);
       // off diagonal entries are a/2, from y = ax + b
       double c = a2/4.0;
       // diagonal coefficient is b, from y = ax + b
       double d0 = 2.0*(mUpper + mLower)/a2;
       for (std::size_t col = 0; col < columns ; ++col)
       {
          // linear index (:,n,k)
          std::size_t nko = Nout*col;
          std::size_t nki = Nin*col;

          // Set aliasing modes to zero
          for (std::size_t n = nDealias; n < Nout; ++n)
          {
             out.data()[nko+n] = 0;
          }

          typename Tout::ScalarType* ptr = in.data() + nki;
          for(std::size_t ip = 0; ip < ps.size(); ip++)
          {
             const std::size_t& p = ps[ip];
             const std::size_t& t = ts[ip];

             for(std::size_t i = 0; i < p; i++)
             {
                // shift rhs and zero last
                std::size_t k = 0;
                if(t == 0)
                {
                   ptr++;
                   for (; k < nDealias-1; ++k, ++ptr)
                   {
                      out.data()[nko+k] = *ptr;
                   }
                   out.data()[nko + k] = 0;

                   ptr = out.data() + nko;
                }
                // multiply by Y, shift and zero last
                else
                {
                   std::size_t s;
                   for(std::size_t j = 1; j <= t; j++)
                   {
                      k = 0;
                      typename Tout::ScalarType vPrev = *ptr;
                      typename Tout::ScalarType vCurr = *ptr;
                      if(j == t)
                      {
                         s = 1;
                         ptr++;
                      }
                      else
                      {
                         s = 0;
                         out.data()[nko+k] = c*((*ptr)*d0 + *(ptr+1)*2.0);
                         k++; ptr++;
                      }

                      for (; k < nDealias-1-s; ++k, ++ptr)
                      {
                         vCurr = *ptr;
                         out.data()[nko+k] = c*(vPrev + vCurr*d0 + *(ptr+1));
                         vPrev = vCurr;
                      }
                      vCurr = *ptr;
                      out.data()[nko+k] = c*(vPrev + vCurr*d0);
                      vPrev = vCurr;
                      k++;
                      out.data()[nko+k] = c*(vPrev);
                      k++;

                      ptr = out.data() + nko;
                   }
                }

                // Compute derivative
                for (; k > 0; --k)
                {
                   double scale = static_cast<double>(4*k)*a2;
                   out.data()[nko+k-1] = out.data()[nko+k+1] + scale*out.data()[nko+k-1];
                }
             }
          }
       }
    }
}

// explicit instantations
template class SpecOp<mods_t, mods_t, spec_id>;
template class SpecOp<mods_t, mods_t, spec_pad>;
template class SpecOp<mods_t, mods_t, spec_d1>;
template class SpecOp<mods_t, mods_t, spec_d2>;
template class SpecOp<mods_t, mods_t, spec_d3>;
template class SpecOp<mods_t, mods_t, spec_d4>;
template class SpecOp<mods_t, mods_t, spec_d1y1>;
template class SpecOp<mods_t, mods_t, spec_d1y2d1>;

} // namespace Cpu
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
