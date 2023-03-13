/**
 * @file FftRtoC.cpp
 * @brief Fftw Real to Complex backend
 */

// External includes
//
#include <fftw3.h>
#include <cassert>
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

template<class AttIn, class AttOut>
FftOp<View<std::complex<double>, AttOut>, View<double, AttIn>>::FftOp()
{
    // FFTW Fixture
    Library::getInstance();
}

template<class AttIn, class AttOut>
FftOp<View<std::complex<double>, AttOut>, View<double, AttIn>>::~FftOp()
{
    // Destroy plan
    if(_plan != nullptr)
    {
        fftw_destroy_plan(static_cast<fftw_plan>(_plan));
        _plan = nullptr;
    }
}

namespace details
{
    fftw_plan setPlanRtoC(const int fwdSize, const int blockSize)
    {
        using fwdType = double;
        using bwdType = fftw_complex;

        // create temporary storage for plan computation
        const int bwdSize = std::floor(fwdSize/2) + 1;
        std::vector<fwdType> fwdTmp(fwdSize*blockSize);
        std::vector<bwdType> bwdTmp(bwdSize*blockSize);

        const int *fftSize = &fwdSize;

        // Create the complex to real plan
        auto fftwPlan = fftw_plan_many_dft_r2c(1, fftSize, blockSize, fwdTmp.data(), NULL,
            1, fwdSize, bwdTmp.data(), NULL, 1, bwdSize, Library::planFlag());
        if(fftwPlan == NULL)
        {
            throw  std::logic_error("FFTW plan failed!");
        }
        return fftwPlan;
    }
}


template<class AttIn, class AttOut>
void FftOp<View<std::complex<double>, AttOut>, View<double, AttIn>>::applyImpl(View<std::complex<double>, AttOut>& mods, const View<double, AttIn>& phys)
{
    assert(std::floor(phys.dims()[0]/2) + 1 == mods.dims()[0]);
    if(_plan == nullptr)
    {
        Profiler::RegionFixture<4> fix("FftOp::initFft");
        int columns = 0;
        if constexpr(std::is_same_v<AttIn, dense2D>)
        {
            assert(phys.dims()[1] == mods.dims()[1]);
            columns = phys.dims()[1];
        }
        else if constexpr(std::is_same_v<AttIn, DCCSC3D>)
        {
            assert(phys.indices()[1].size() == mods.indices()[1].size());
            columns = phys.indices()[1].size();
        }
        else
        {
            throw std::logic_error("Not implemented yet.");
        }
        _plan = details::setPlanRtoC(phys.dims()[0], columns);
    }
    Profiler::RegionFixture<4> fix("FftOp::applyFft");
    fftw_execute_dft_r2c(static_cast<fftw_plan>(_plan),
        const_cast<double* >(phys.data()),
        reinterpret_cast<fftw_complex* >(mods.data()));
}

// Explicit instantiations
template class FftOp<CmodsDense2D_t, RphysDense2D_t>;
template class FftOp<CmodsDCCSC3D_t, RphysDCCSC3D_t>;


} // namespace Fftw
} // namespace Fft
} // namespace QuICC
