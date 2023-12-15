/**
 * @file FftCtoC.cpp
 * @brief CuFft Complex to Complex backend
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
FftOp<View::View<std::complex<double>, AttOut>, View::View<std::complex<double>, AttIn>>::~FftOp()
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
    cufftHandle* setPlanCtoC(const int fwdSize, const int blockSize)
    {
        // Create the complex to real plan
        cufftHandle* plan = new cufftHandle;
        constexpr int rank = 1;
        constexpr int istride = 1;
        constexpr int ostride = 1;
        const int bwdSize = fwdSize;
        auto err = cufftPlanMany(plan, rank, const_cast<int*>(&fwdSize), NULL, istride,
            bwdSize, NULL, ostride, fwdSize, CUFFT_Z2Z, blockSize);
        if(err != CUFFT_SUCCESS)
        {
            throw  std::logic_error("CuFFT plan failed! error code: "+std::to_string(err) );
        }
        return plan;
    }
}

template<class AttIn, class AttOut>
void FftOp<View::View<std::complex<double>, AttOut>, View::View<std::complex<double>, AttIn>>::applyImpl(View::View<std::complex<double>, AttOut>& phys, const View::View<std::complex<double>, AttIn>& mods)
{
    using namespace QuICC::View;
    if(_plan == nullptr)
    {
        Profiler::RegionFixture<5> fix("FftOp::initFft");
        int columns = 0;
        if constexpr(std::is_same_v<AttIn, dense2D>)
        {
            assert(phys.dims()[0] == mods.dims()[0]);
            assert(phys.dims()[1] == mods.dims()[1]);
            columns = phys.dims()[1];
        }
        else if constexpr(std::is_same_v<AttIn, DCCSC3DInOrder>)
        {
            assert(phys.dims()[0] == mods.lds());
            assert(phys.indices()[1].size() == mods.indices()[1].size());
            columns = phys.indices()[1].size();
        }
        else
        {
            throw std::logic_error("Not implemented yet.");
        }
        _plan = details::setPlanCtoC(phys.dims()[0], columns);
    }
    Profiler::RegionFixture<5> fix("FftOp::applyFft");
    auto err = cufftExecZ2Z(*static_cast<cufftHandle*>(_plan),
        reinterpret_cast<cufftDoubleComplex*>(mods.data()),
        reinterpret_cast<cufftDoubleComplex*>(phys.data()), CUFFT_INVERSE);
    if(err != CUFFT_SUCCESS)
    {
        throw  std::logic_error("CuFFT execute failed! error code: "+std::to_string(err) );
    }
}

// Explicit instantiations
template class FftOp<CphysDense2D_t, CmodsDense2D_t>;
template class FftOp<CphysDCCSC3DInOrder_t, CmodsDCCSC3DInOrder_t>;


} // namespace CuFft
} // namespace Fft
} // namespace QuICC
