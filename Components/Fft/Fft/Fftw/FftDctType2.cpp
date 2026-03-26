/**
 * @file FftDctType2.cpp
 * @brief Fftw DCT Type 2 backend
 */

// External includes
//
#include <cassert>
#include <fftw3.h>
#include <stdexcept>
#include <type_traits>

// Project includes
//
#include "Fft.hpp"
#include "Fft/FftTypes.hpp"
#include "Library.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Fft {
namespace Fftw {

namespace details {
fftw_plan setPlanDctType2(const int fwdSize, const int blockSize, const int mult);
} // namespace details

//
// DCT of real data
//

template <class AttIn, class AttOut>
FftOp<View::View<double, AttOut>,
   View::View<double, AttIn>, dct_type2_t>::FftOp()
{
   // FFTW Fixture
   Library::getInstance();
}

template <class AttIn, class AttOut>
FftOp<View::View<double, AttOut>,
   View::View<double, AttIn>, dct_type2_t>::~FftOp()
{
   // Destroy plan
   if (_plan != nullptr)
   {
      fftw_destroy_plan(static_cast<fftw_plan>(_plan));
      _plan = nullptr;
   }
}


template <class AttIn, class AttOut>
void FftOp<View::View<double, AttOut>,
   View::View<double, AttIn>, dct_type2_t>::applyImpl(View::View<double,
                                            AttOut>& mods,
   const View::View<double, AttIn>& phys)
{
   using namespace QuICC::View;
   if (_plan == nullptr)
   {
      Profiler::RegionFixture<5> fix("Fftw::FftOp::initFft-DctType2");
      int columns = 0;
      if constexpr (std::is_same_v<AttIn, dense2D>)
      {
         assert(phys.dims()[0] == mods.dims()[0]);
         assert(phys.dims()[1] == mods.dims()[1]);
         columns = phys.dims()[1];
      }
      else if constexpr (std::is_same_v<AttIn, DCCSC3D>)
      {
         assert(phys.dims()[0] == mods.lds());
         assert(phys.indices()[1].size() == mods.indices()[1].size());
         columns = phys.indices()[1].size();
      }
      else
      {
         throw std::logic_error("Not implemented yet.");
      }
      _plan = details::setPlanDctType2(phys.dims()[0], columns, 1);
   }
   Profiler::RegionFixture<5> fix("Fftw::FftOp::applyFft-DctType2");
   fftw_execute_r2r(static_cast<fftw_plan>(_plan),
      const_cast<double*>(phys.data()),
      mods.data());
}

// Explicit instantiations
template class FftOp<RphysDense2D_t, RmodsDense2D_t, dct_type2_t>;
template class FftOp<RphysDCCSC3D_t, RmodsDCCSC3D_t, dct_type2_t>;

//
// DCT of complex data
//

template <class AttIn, class AttOut>
FftOp<View::View<std::complex<double>, AttOut>,
   View::View<std::complex<double>, AttIn>, dct_type2_t>::FftOp()
{
   // FFTW Fixture
   Library::getInstance();
}

template <class AttIn, class AttOut>
FftOp<View::View<std::complex<double>, AttOut>,
   View::View<std::complex<double>, AttIn>, dct_type2_t>::~FftOp()
{
   // Destroy plan
   if (_plan != nullptr)
   {
      fftw_destroy_plan(static_cast<fftw_plan>(_plan));
      _plan = nullptr;
   }
}


template <class AttIn, class AttOut>
void FftOp<View::View<std::complex<double>, AttOut>,
   View::View<std::complex<double>, AttIn>, dct_type2_t>::applyImpl(View::View<std::complex<double>,
                                            AttOut>& mods,
   const View::View<std::complex<double>, AttIn>& phys)
{
   using namespace QuICC::View;
   if (_plan == nullptr)
   {
      Profiler::RegionFixture<5> fix("Fftw::FftOp::initFft-DctType2");
      int columns = 0;
      if constexpr (std::is_same_v<AttIn, dense2D>)
      {
         assert(phys.dims()[0] == mods.dims()[0]);
         assert(phys.dims()[1] == mods.dims()[1]);
         columns = phys.dims()[1];
      }
      else if constexpr (std::is_same_v<AttIn, DCCSC3D>)
      {
         assert(phys.dims()[0] == mods.lds());
         assert(phys.indices()[1].size() == mods.indices()[1].size());
         columns = phys.indices()[1].size();
      }
      else
      {
         throw std::logic_error("Not implemented yet.");
      }
      _plan = details::setPlanDctType2(phys.dims()[0], columns, 2);
   }
   Profiler::RegionFixture<5> fix("Fftw::FftOp::applyFft-DctType2");
   fftw_execute_r2r(static_cast<fftw_plan>(_plan),
      const_cast<double*>(reinterpret_cast<double*>(phys.data())),
      reinterpret_cast<double*>(mods.data()));
   fftw_execute_r2r(static_cast<fftw_plan>(_plan),
      const_cast<double*>(reinterpret_cast<double*>(phys.data()))+1,
      reinterpret_cast<double*>(mods.data())+1);
}

// Explicit instantiations
template class FftOp<CphysDense2D_t, CmodsDense2D_t, dct_type2_t>;
template class FftOp<CphysDCCSC3D_t, CmodsDCCSC3D_t, dct_type2_t>;


namespace details {
fftw_plan setPlanDctType2(const int fwdSize, const int blockSize, const int mult)
{
   using fwdType = double;
   using bwdType = double;

   // create temporary storage for plan computation
   const int bwdSize = fwdSize;
   std::vector<fwdType> fwdTmp(mult*fwdSize * blockSize);
   std::vector<bwdType> bwdTmp(mult*bwdSize * blockSize);

   const int* fftSize = &fwdSize;

   // Create the real to real type II plan
   const fftw_r2r_kind fftKind[] = {FFTW_REDFT10};
   auto fftwPlan = fftw_plan_many_r2r(1, fftSize, blockSize, bwdTmp.data(),
      NULL, mult, mult*bwdSize, fwdTmp.data(), NULL, mult, mult*fwdSize, fftKind, Library::planFlag());
   if (fftwPlan == NULL)
   {
      throw std::logic_error("FFTW plan failed!");
   }
   return fftwPlan;
}
} // namespace details

} // namespace Fftw
} // namespace Fft
} // namespace QuICC
