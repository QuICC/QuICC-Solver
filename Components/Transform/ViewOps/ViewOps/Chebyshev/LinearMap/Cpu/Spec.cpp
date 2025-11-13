#include <complex>
#include <memory>

#include "Spec.hpp"
#include "Profiler/Interface.hpp"
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
SpecOp<Tout, Tin, Operation, Treatment>::SpecOp(const double lower,
   const double upper) :
    mLower(lower), mUpper(upper), mScale(1.0){};

template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
SpecOp<Tout, Tin, Operation, Treatment>::SpecOp(const double lower,
   const double upper, ScaleType scale) :
    mLower(lower), mUpper(upper), mScale(scale){};

template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
void SpecOp<Tout, Tin, Operation, Treatment>::applyImpl(Tout& out,
   const Tin& in, const ScaleType fftScaling)
{
   Profiler::RegionFixture<5> fix("SpecOp::applyImpl");

#ifdef QUICC_HAS_CUDA_BACKEND
   assert(!QuICC::Cuda::isDeviceMemory(out.data()));
#endif

   assert(out.dims()[1] == in.dims()[1]);
   assert(out.dims()[2] == in.dims()[2]);

   if constexpr (std::is_same_v<Operation, spec_id>)
   {
      // if the spec is identity and in place and there are no modes
      // to be zeroed then it is a noop
      if (out.data() == in.data() && out.dims()[0] == out.lds())
      {
         return;
      }
   }

   std::size_t nDealias;
   if constexpr (Treatment & ndealias_out)
   {
      assert(out.size() <= in.size());
      assert(out.dims()[0] <= in.dims()[0]);

      // if the spec is identity and in place and there are no modes
      // to be zeroed then it is a noop
      if (out.data() == in.data() && out.dims()[0] == out.lds())
      {
         return;
      }

      nDealias = out.dims()[0];
   }
   else if constexpr (Treatment & ndealias_in)
   {
      assert(out.size() >= in.size());
      assert(out.dims()[0] >= in.dims()[0]);

      nDealias = in.dims()[0];
   }
   else if constexpr (Treatment == none_t)
   {
      assert(out.size() >= in.size());
      assert(out.dims()[0] >= in.dims()[0]);

      nDealias = out.lds();
   }
   else
   {
      throw std::logic_error("Unknown Treatment parameter");
   }

   std::size_t Nout = out.lds();
   std::size_t Nin = in.lds();

   ScaleType c = fftScaling;

   // Column major
   // Get total number of columns to loop over
   auto indices = in.indices()[1];
   auto columns = indices.size();

   std::size_t col0 = 0;
   std::size_t mean_cols = 0;
   if constexpr (Treatment & (zero_l0 | mean_op))
   {
      mean_cols = out.pointers()[1][1] - out.pointers()[1][0];
   }

   if constexpr (Treatment & zero_l0)
   {
      for (std::size_t col = 0; col < mean_cols; ++col)
      {
         std::size_t nko = Nout * col;
         for (std::size_t n = 0; n < nDealias; ++n)
         {
            out.data()[nko + n] = 0;
         }
      }
      col0 = mean_cols;
   }

   if constexpr (Operation::p.size() == 1 && Operation::p[0] == 0)
   {
      const std::size_t& t = Operation::t[0];
      for (std::size_t col = col0; col < columns; ++col)
      {
         // linear index (:,n,k)
         std::size_t nko = Nout * col;
         std::size_t nki = Nin * col;
         std::size_t n = 0;

         if (t == 0)
         {
            if constexpr(std::is_same_v<Operation, spec_int>)
            {
               ScaleType cc = c * 2.0 * (this->mUpper - this->mLower);
               for (; n < nDealias; ++n)
               {
                  ScaleType dn = static_cast<ScaleType>(n);
                  out.data()[nko + n] = in.data()[nki + n] * cc /(1 - dn*dn);
                  n++;
                  out.data()[nko + n] *= 0;
               }
               out.data()[nko] *= 0.5;
            }
            else
            {
               for (; n < nDealias; ++n)
               {
                  out.data()[nko + n] = in.data()[nki + n] * c;
               }
            }
         }
         else
         {
            if constexpr (Treatment & ndealias_out)
            {
               n = nDealias + 1;
            }
            else if constexpr (Treatment & ndealias_in)
            {
               n = in.dims()[0];
            }
            typename Tout::ScalarType* inPtr = in.data() + nki;
            typename Tout::ScalarType* outPtr = out.data() + nko;
            for (std::size_t j = 1; j <= t; j++)
            {
               n = this->multiplyByY(outPtr, inPtr, Nout, n, c);
               inPtr = out.data() + nko;
            }
         }

         if constexpr (Treatment & zero_pad)
         {
            for (; n < Nout; ++n)
            {
               out.data()[nko + n] = 0;
            }
         }
      }
   }
   else
   {
      const auto& ps = Operation::p;
      const auto& ts = Operation::t;

      for (std::size_t col = col0; col < columns; ++col)
      {
         // linear index (:,n,k)
         std::size_t nko = Nout * col;
         std::size_t nki = Nin * col;

         typename Tout::ScalarType* inPtr = in.data() + nki;
         typename Tout::ScalarType* outPtr = out.data() + nko;

         // Set aliasing modes to zero
         if constexpr (Treatment & zero_pad)
         {
            for (std::size_t n = nDealias; n < Nout; ++n)
            {
               outPtr[n] = 0;
            }
         }

         for (std::size_t ip = 0; ip < ps.size(); ip++)
         {
            const std::size_t& p = ps[ip];
            const std::size_t& t = ts[ip];

            for (std::size_t i = 0; i < p; i++)
            {
               // shift rhs and zero last
               std::size_t k = 0;
               if (t == 0)
               {
                  inPtr++;
                  for (; k < nDealias - 1; ++k, ++inPtr)
                  {
                     outPtr[k] = *inPtr;
                  }
                  outPtr[k] = 0;

                  inPtr = out.data() + nko;
               }
               // multiply by Y, shift and zero last
               else
               {
                  k = nDealias;
                  std::size_t s = 0;
                  for (std::size_t j = 1; j <= t; j++)
                  {
                     if (j == t)
                     {
                        s = 1;
                     }
                     k = this->multiplyByY(outPtr, inPtr, Nout, k - s, c, s);
                     inPtr = out.data() + nko;
                  }
               }

               // Compute derivative
               differentiate(outPtr, k + 1);
            }
         }
      }
   }

   if constexpr(Treatment & mean_op)
   {
      if(!std::is_same_v<typename Operation::SparseMeanOpType, void>)
      {
         if(!this->mSparseMeanOp)
         {
            typename Operation::SparseMeanOpType op(nDealias, nDealias, this->mLower, this->mUpper);
            this->mSparseMeanOp = std::make_unique<SparseMatrix>(op.mat());
         }

         Eigen::Map<MatrixZ, Eigen::Unaligned, Eigen::Stride<::Eigen::Dynamic,1> > outMap(out.data(), nDealias, mean_cols, Eigen::Stride<::Eigen::Dynamic,1>(Nout, 1));

         outMap = (*this->mSparseMeanOp) * outMap;
      }
   }

   if constexpr(!std::is_same_v<typename Operation::SparseOpType, void>)
   {
      if(!this->mSparseOp)
      {
         typename Operation::SparseOpType op(nDealias, nDealias, this->mLower, this->mUpper);
         this->mSparseOp = std::make_unique<SparseMatrix>(op.mat());
      }

      Eigen::Map<MatrixZ, Eigen::Unaligned, Eigen::Stride<::Eigen::Dynamic,1> > outMap(out.data() + Nout*mean_cols, nDealias, columns - mean_cols, Eigen::Stride<::Eigen::Dynamic,1>(Nout, 1));

      outMap = (*this->mSparseOp) * outMap;
   }
}

template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
std::size_t SpecOp<Tout, Tin, Operation, Treatment>::multiplyByY(
   typename Tout::ScalarType* const out, typename Tout::ScalarType* in,
   const std::size_t Nout, const std::size_t Nin, const double c,
   const std::size_t shiftIn)
{
   assert(Nout >= Nin - 2);

   // 2*a, from y=ax + b
   double a2 = (mUpper - mLower);
   // off diagonal entries are a/2, from y = ax + b
   double d1 = c * a2 / 4.0;
   // diagonal coefficient is b, from y = ax + b
   double d0 = 2.0 * (mUpper + mLower) / a2;

   std::size_t k;
   const std::size_t Nk = Nin - 1;
   typename Tout::ScalarType vPrev = in[0];
   typename Tout::ScalarType vCurr;

   if (shiftIn > 0)
   {
      k = 0;
      in += shiftIn;
   }
   else
   {
      out[0] = d1 * (in[0] * d0 + in[1] * 2.0);
      k = 1;
      in++;
   }

   for (; k < Nk; ++k, ++in)
   {
      vCurr = in[0];
      out[k] = d1 * (vPrev + in[0] * d0 + in[1]);
      vPrev = vCurr;
   }
   if (k < Nout)
   {
      vCurr = in[0];
      out[k] = d1 * (vPrev + in[0] * d0);
      vPrev = vCurr;
      k++;
   }
   if (k < Nout)
   {
      out[k] = d1 * (vPrev);
      k++;
   }

   return k;
}

template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
void SpecOp<Tout, Tin, Operation, Treatment>::differentiate(
   typename Tout::ScalarType* const out, const std::size_t Nout)
{
   // Compute derivative
   assert(out[Nout] == 0.0);
   if constexpr(std::is_same_v<std::complex<double>,typename Tout::ScalarType>)
   {
      // 2*a, from y=ax + b
      double c = 2.0*(mUpper - mLower);

      double* dout = reinterpret_cast<double *>(out);
      for (std::size_t k = 2*(Nout - 1); k > 0; k -= 2)
      {
         double scale = static_cast<double>(k) * c;
         dout[k - 2] = dout[k + 2] + scale * dout[k - 2];
         dout[k - 1] = dout[k + 3] + scale * dout[k - 1];
      }
   }
   else
   {
      // 2*a, from y=ax + b
      double c = 4.0*(mUpper - mLower);

      for (std::size_t k = Nout - 1; k > 0; --k)
      {
         double scale = static_cast<double>(k) * c;
         out[k - 1] = out[k + 1] + scale * out[k - 1];
      }
   }
}

// explicit instantations
template class SpecOp<mods_t, mods_t, spec_id, ndealias_out>;
template class SpecOp<mods_t, mods_t, spec_id, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_y1, ndealias_out>;
template class SpecOp<mods_t, mods_t, spec_y1, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_i2, ndealias_out>;
template class SpecOp<mods_t, mods_t, spec_i2, ndealias_out | mean_op>;
template class SpecOp<mods_t, mods_t, spec_i2d1, ndealias_out>;
template class SpecOp<mods_t, mods_t, spec_i2d1, ndealias_out | mean_op>;
template class SpecOp<mods_t, mods_t, spec_i2y1d1y1, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_i2y1, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_i2y2d1y1, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_i2y2, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_i4, ndealias_out>;
template class SpecOp<mods_t, mods_t, spec_i4y3d1y1, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_i4y3, ndealias_out | zero_l0>;
template class SpecOp<mods_t, mods_t, spec_i4d1, ndealias_out>;
template class SpecOp<mods_t, mods_t, spec_i4d1, ndealias_out | mean_op>;
template class SpecOp<mods_t, mods_t, spec_id, ndealias_in | zero_pad>;
template class SpecOp<mods_t, mods_t, spec_y1, ndealias_in | zero_pad>;
template class SpecOp<mods_t, mods_t, spec_d1, ndealias_in | zero_pad>;
template class SpecOp<mods_t, mods_t, spec_d2, ndealias_in | zero_pad>;
template class SpecOp<mods_t, mods_t, spec_d3, ndealias_in | zero_pad>;
template class SpecOp<mods_t, mods_t, spec_d4, ndealias_in | zero_pad>;
template class SpecOp<mods_t, mods_t, spec_d1y1, ndealias_in | zero_pad>;
template class SpecOp<mods_t, mods_t, spec_d1y2d1, ndealias_in | zero_pad>;
template class SpecOp<power_t, power_t, spec_int, none_t>;

} // namespace Cpu
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
