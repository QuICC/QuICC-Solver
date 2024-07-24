/**
 * @file FftCtoR.cpp
 * @brief CuFft Complex to Real backend
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
FftOp<View::View<double, AttOut>, View::View<std::complex<double>, AttIn>>::~FftOp()
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
    cufftHandle* setPlanCtoR(const int fwdSize, const int blockSize)
    {
        // Create the complex to real plan
        cufftHandle* plan = new cufftHandle;
        constexpr int rank = 1;
        constexpr int istride = 1;
        constexpr int ostride = 1;
        const int idist = std::floor(fwdSize/2) + 1;
        const int odist = fwdSize;
        auto err = cufftPlanMany(plan, rank, const_cast<int*>(&fwdSize),
            NULL, istride, idist,
            NULL, ostride, odist,
            CUFFT_Z2D, blockSize);
        if(err != CUFFT_SUCCESS)
        {
            throw  std::logic_error("CuFFT plan failed! error code: "+std::to_string(err) );
        }
        return plan;
    }
}

template<class AttIn, class AttOut>
void FftOp<View::View<double, AttOut>, View::View<std::complex<double>, AttIn>>::applyImpl(View::View<double, AttOut>& phys, const View::View<std::complex<double>, AttIn>& mods)
{
    using namespace QuICC::View;
    if(_plan == nullptr)
    {
        Profiler::RegionFixture<5> fix("FftOp::initFft");
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
    Profiler::RegionFixture<5> fix("FftOp::applyFft");
    auto err = cufftExecZ2D(*static_cast<cufftHandle*>(_plan),
        reinterpret_cast<cufftDoubleComplex*>(mods.data()),
        phys.data());
    if(err != CUFFT_SUCCESS)
    {
        throw  std::logic_error("CuFFT execute failed! error code: "+std::to_string(err) );
    }
}

// Explicit instantiations
template class FftOp<RphysDense2D_t, CmodsDense2D_t>;
template class FftOp<RphysDCCSC3D_t, CmodsDCCSC3D_t>;


} // namespace CuFft
} // namespace Fft
} // namespace QuICC
