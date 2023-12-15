/**
 * @file FftCtoR.cpp
 * @brief Fftw Complex to Real backend
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
FftOp<View::View<double, AttOut>, View::View<std::complex<double>, AttIn>>::FftOp()
{
    // FFTW Fixture
    Library::getInstance();
}

template<class AttIn, class AttOut>
FftOp<View::View<double, AttOut>, View::View<std::complex<double>, AttIn>>::~FftOp()
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
    fftw_plan setPlanCtoR(const int fwdSize, const int blockSize)
    {
        using fwdType = double;
        using bwdType = fftw_complex;

        // create temporary storage for plan computation
        const int bwdSize = std::floor(fwdSize/2) + 1;
        std::vector<fwdType> fwdTmp(fwdSize*blockSize);
        std::vector<bwdType> bwdTmp(bwdSize*blockSize);

        const int *fftSize = &fwdSize;

        // Create the complex to real plan
        auto fftwPlan = fftw_plan_many_dft_c2r(1, fftSize, blockSize, bwdTmp.data(), NULL,
            1, bwdSize, fwdTmp.data(), NULL, 1, fwdSize, Library::planFlag());
        if(fftwPlan == NULL)
        {
            throw  std::logic_error("FFTW plan failed!");
        }
        return fftwPlan;
    }
}


template<class AttIn, class AttOut>
void FftOp<View::View<double, AttOut>, View::View<std::complex<double>, AttIn>>::applyImpl(View::View<double, AttOut>& phys, const View::View<std::complex<double>, AttIn>& mods)
{
    using namespace QuICC::View;
    if(_plan == nullptr)
    {
        Profiler::RegionFixture<5> fix("Fftw::FftOp::initFft-CtoR");
        int columns = 0;
        if constexpr(std::is_same_v<AttIn, dense2D>)
        {
            assert(std::floor(phys.dims()[0]/2) + 1 == mods.dims()[0]);
            assert(phys.dims()[1] == mods.dims()[1]);
            columns = phys.dims()[1];
        }
        else if constexpr(std::is_same_v<AttIn, DCCSC3D>)
        {
            assert(std::floor(phys.dims()[0]/2) + 1 == mods.lds());
            assert(phys.indices()[1].size() == mods.indices()[1].size());
            columns = phys.indices()[1].size();
        }
        else
        {
            throw std::logic_error("Not implemented yet.");
        }
        _plan = details::setPlanCtoR(phys.dims()[0], columns);
    }
    Profiler::RegionFixture<5> fix("Fftw::FftOp::applyFft-CtoR");
    fftw_execute_dft_c2r(static_cast<fftw_plan>(_plan),
        reinterpret_cast<fftw_complex* >(
        const_cast<std::complex<double>*>(mods.data())), phys.data());
}

// Explicit instantiations
template class FftOp<RphysDense2D_t, CmodsDense2D_t>;
template class FftOp<RphysDCCSC3D_t, CmodsDCCSC3D_t>;

} // namespace Fftw
} // namespace Fft
} // namespace QuICC
