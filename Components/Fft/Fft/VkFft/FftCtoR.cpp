/**
 * @file FftCtoR.cpp
 * @brief VkFft Complex to Real backend
 */

// External includes
//
#include <vkFFT.h>

// Project includes
//
#include "Fft.hpp"
#include "Fft/FftTypes.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Fft {
namespace VkFft {

template <class AttIn, class AttOut>
FftOp<View::View<double, AttOut>,
   View::View<std::complex<double>, AttIn>>::~FftOp()
{
   // Destroy plan
   if (_plan != nullptr)
   {
      free(((VkFFTApplication*)_plan)->configuration.device);
      deleteVkFFT((VkFFTApplication*)_plan);
      delete static_cast<VkFFTApplication*>(_plan);
      _plan = nullptr;
   }
}

namespace details {
VkFFTApplication* setPlanCtoR(const int fwdSize, const int blockSize)
{
   // Create the complex to real plan
   VkFFTApplication* plan = new VkFFTApplication();

   VkFFTConfiguration configuration = VKFFT_ZERO_INIT;
   configuration.FFTdim = 1;

   FILE* kernelCache;
   uint64_t str_len;
   char fname[500];
#ifdef VKFFT_USE_DOUBLEDOUBLE_FP128
   sprintf(fname, "%s/VkFFT_%d_%d_%d_c2r_qddp", VKFFT_KERNELS_DIR,
      VkFFTGetVersion(), fwdSize, blockSize);
#else
   sprintf(fname, "%s/VkFFT_%d_%d_%d_c2r_dp", VKFFT_KERNELS_DIR,
      VkFFTGetVersion(), fwdSize, blockSize);
#endif
   kernelCache = fopen(fname, "rb");
   if (kernelCache != 0)
   {
      fseek(kernelCache, 0, SEEK_END);
      str_len = ftell(kernelCache);
      fseek(kernelCache, 0, SEEK_SET);
      if (str_len < 10)
      {
         configuration.saveApplicationToString = 0;
         fclose(kernelCache);
      }
      else
      {
         configuration.loadApplicationFromString = 1;
         configuration.loadApplicationString = malloc(str_len);
         fread(configuration.loadApplicationString, str_len, 1, kernelCache);
         fclose(kernelCache);
      }
   }
   else
   {
      configuration.saveApplicationToString = 1;
   }

   configuration.size[0] = fwdSize;
   configuration.numberBatches = blockSize;
#ifdef VKFFT_USE_DOUBLEDOUBLE_FP128
   configuration.quadDoubleDoublePrecisionDoubleMemory = 1;
#else
   configuration.doublePrecision = 1;
#endif
   configuration.isOutputFormatted = 1;
   configuration.performR2C = 1;
   configuration.makeInversePlanOnly = 1;

   CUresult err = CUDA_SUCCESS;
   err = cuInit(0);
   if (err != CUDA_SUCCESS)
   {
      throw std::logic_error("CUDA failed! error code: " + std::to_string(err));
   }
   configuration.device = (CUdevice*)calloc(1, sizeof(CUdevice));
   err = cuDeviceGet(configuration.device,
      0); // need to fix id for multi-GPU/node systems
   if (err != CUDA_SUCCESS)
   {
      throw std::logic_error("CUDA failed! error code: " + std::to_string(err));
   }
   VkFFTResult resFFT = initializeVkFFT(plan, configuration);

   if (configuration.loadApplicationFromString)
   {
      free(configuration.loadApplicationString);
   }

   if (configuration.saveApplicationToString)
   {
      kernelCache = fopen(fname, "wb");
      fwrite(plan->saveApplicationString, plan->applicationStringSize, 1,
         kernelCache);
      fclose(kernelCache);
   }

   if (resFFT != VKFFT_SUCCESS)
   {
      std::string str(getVkFFTErrorString(resFFT));
      throw std::logic_error("VkFFT plan failed! error code: " + str);
   }
   return plan;
}
} // namespace details

template <class AttIn, class AttOut>
void FftOp<View::View<double, AttOut>,
   View::View<std::complex<double>, AttIn>>::applyImpl(View::View<double,
                                                          AttOut>& phys,
   const View::View<std::complex<double>, AttIn>& mods)
{
   using namespace QuICC::View;
   if (_plan == nullptr)
   {
      Profiler::RegionFixture<5> fix("FftOp::initFft");
      int columns = 0;
      if constexpr (std::is_same_v<AttIn, dense2D>)
      {
         assert(std::floor(phys.dims()[0] / 2) + 1 == mods.dims()[0]);
         assert(phys.dims()[1] == mods.dims()[1]);
         columns = phys.dims()[1];
      }
      else if constexpr (std::is_same_v<AttIn, DCCSC3D>)
      {
         assert(std::floor(phys.dims()[0] / 2) + 1 == mods.lds());
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

   VkFFTLaunchParams launchParams = VKFFT_ZERO_INIT;
   void* inputData = mods.data();
   void* outputData = phys.data();
   launchParams.buffer = &inputData;
   launchParams.outputBuffer = &outputData;
   VkFFTResult resFFT = VkFFTAppend((VkFFTApplication*)_plan, 1, &launchParams);
   if (resFFT != VKFFT_SUCCESS)
   {
      std::string str(getVkFFTErrorString(resFFT));
      throw std::logic_error("VkFFT execute failed! error code: " + str);
   }
}

// Explicit instantiations
template class FftOp<RphysDense2D_t, CmodsDense2D_t>;
template class FftOp<RphysDCCSC3D_t, CmodsDCCSC3D_t>;


} // namespace VkFft
} // namespace Fft
} // namespace QuICC
