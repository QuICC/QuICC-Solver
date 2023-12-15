/**
 * @file FftRtoC.cpp
 * @brief CuFft Real to Complex backend
 */

// External includes
//
#include <iostream>
#include <cufft.h>

// Project includes
//
#include "Fft.hpp"
#include "Fft/FftTypes.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Fft {
namespace CuFft {

template<class AttIn, class AttOut>
FftOp<View::View<std::complex<double>, AttOut>, View::View<double, AttIn>>::~FftOp()
{
    // Destroy plan
    if(_plan != nullptr)
    {
        auto err = cufftDestroy(*static_cast<cufftHandle*>(_plan));
        if(err != CUFFT_SUCCESS)
        {
            std::cerr << "CuFFT destroy plan failed! error code: "
                +std::to_string(err)+"\n";
        }
        delete static_cast<cufftHandle*>(_plan);
        _plan = nullptr;
    }
}

namespace details
{
    cufftHandle* setPlanRtoC(const int fwdSize, const int blockSize)
    {
        // Create the real to complex plan
        cufftHandle* plan = new cufftHandle;
        constexpr int rank = 1;
        constexpr int istride = 1;
        constexpr int ostride = 1;
        const int idist = fwdSize;
        const int odist = std::floor(fwdSize/2) + 1;
        auto err = cufftPlanMany(plan, rank, const_cast<int*>(&fwdSize),
            NULL, istride, idist,
            NULL, ostride, odist,
            CUFFT_D2Z, blockSize);
        if(err != CUFFT_SUCCESS)
        {
            throw  std::logic_error("CuFFT plan failed! error code: "+std::to_string(err) );
        }
        return plan;
    }
}

template<class AttIn, class AttOut>
void FftOp<View::View<std::complex<double>, AttOut>, View::View<double, AttIn>>::applyImpl(View::View<std::complex<double>, AttOut>& mods, const View::View<double, AttIn>& phys)
{
    assert(std::floor(phys.dims()[0]/2) + 1 == mods.dims()[0]);
    using namespace QuICC::View;
    if(_plan == nullptr)
    {
        Profiler::RegionFixture<5> fix("FftOp::initFft");
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
    Profiler::RegionFixture<5> fix("FftOp::applyFft");
    auto err = cufftExecD2Z(*static_cast<cufftHandle*>(_plan),
        phys.data(),
        reinterpret_cast<cufftDoubleComplex*>(mods.data()));
    if(err != CUFFT_SUCCESS)
    {
        throw  std::logic_error("CuFFT execute failed! error code: "+std::to_string(err) );
    }
}

// Explicit instantiations
template class FftOp<CmodsDense2D_t, RphysDense2D_t>;
template class FftOp<CmodsDCCSC3D_t, RphysDCCSC3D_t>;


} // namespace CuFft
} // namespace Fft
} // namespace QuICC
