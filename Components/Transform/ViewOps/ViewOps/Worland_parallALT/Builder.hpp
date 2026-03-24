/**
 * @file Builder.hpp
 * @brief Generic Worland parallALT operator builder
 */
#pragma once

// External includes
//
#include <memory>
#include <cstdint>

// Project includes
//
#include "Operator/Binary.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Worland_parallALT/TypeTraits.hpp"

#include <cuComplex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <nvrtc.h>
#include "parallALT.hpp"

namespace QuICC {
namespace Transform {
namespace Worland_parallALT {
using namespace QuICC::Operator;
using type = QuICC::Memory::Cuda::Malloc;
/// @brief Derived classes implement the differentiation in modal space
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Order of differentiation
/// @tparam Direction Fft direction tag
/// @tparam Treatment special treatment mask, typically of mode zero or dealiasing
template<class Tout, class Tin>
class ParallaltOp : public UnaryBaseOp<ParallaltOp<Tout, Tin>, Tout, Tin> {
public:
    /// @brief Default constructor
   ParallaltOp(std::shared_ptr<Memory::memory_resource> mem) : _mem(mem) {
   };
    /// @brief dtor
    ~ParallaltOp()
    {
        if (appContainer.config.Ntheta != 0)
            deleteParallALT(&VkGPU, &appContainer);

        /* if (temp_buffer != 0)
        {
           cudaFree(temp_buffer);
           temp_buffer = 0;
        };
        if (temp_buffer2 != 0)
        {
           cudaFree(temp_buffer2);
           temp_buffer2 = 0;
        };*/
    };

    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    void initImpl(Tout& out, const Tin& in){
        std::uint32_t Ntheta = (transformDirection) ? out.dims()[0] : in.dims()[0];//igrid.size();
        
        //std::uint32_t nLayers = static_cast<std::uint32_t>(this->mspSetup->slowSize());

        ///\todo this should be the full matrix size
        //std::uint32_t M = out.pointers()->size() - 1;

        PfSolve::PfSolveResult resPfSolve = PfSolve::PFSOLVE_SUCCESS;
        config = {};
	    config.Ntheta = Ntheta;//this->mspSetup->bwdSize();
        config.radialTransform = transformType;
        config.useGraphs = 0;
        config.initializeGraph = 0;
        config.useUberKernel = 1;
	    config.testMerge = 0;
	    config.testAccuracy = 1;
	    config.doALTOnly = 1;
	    config.useMatMulConnection = 1;
	    config.use_tc = 2;
	    config.mergeType = 1;
	    config.numMergedIterMatMul = 8;
	    config.numMergedIterMax = 1;
	    config.numMergedIterMin = 1;
	    config.disableCaching = 1;
	    config.fixAccuracy = 1;
	    config.numRadialBatches = 1;
	    config.profile_iter = 1;
	    config.profile_iter_combined = 1;
        config.specifyBuffersAtLaunch = 1;
        config.convertPackedStrided = convertPackedStrided;
	    config.WMMA_M = 8;
	    config.WMMA_N = 8;
	    config.WMMA_K = 4;
	    appContainer = {};
        config.projector = transformDirection;
        
         std::uint32_t* temp_pointers =
           (std::uint32_t*)calloc(out.pointers()[1].size(), sizeof(std::uint32_t));
        cudaErrChk(cudaMemcpyAsync(temp_pointers, out.pointers()[1].data(),
        out.pointers()[1].size() * sizeof(std::uint32_t), cudaMemcpyDeviceToHost));
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
           int numRHS = temp_pointers[i] - temp_pointers[i-1];
           if ((((i - 1) % 2) == 0) && (numRHS!=0))
            {
                config.num_m_even++;
            }
            if ((((i-1)%2) == 1 ) && (numRHS!=0))
            {
                config.num_m_odd++;
            }
        }

        int start = 0;
        int iter = 0;
        int* m_even = (int*)calloc(config.num_m_even, sizeof(int));
        int* m_even_endBatch = (int*)calloc(config.num_m_even, sizeof(int));
        int* m_odd = (int*)calloc(config.num_m_odd, sizeof(int));
        int* m_odd_endBatch = (int*)calloc(config.num_m_odd, sizeof(int));

        int* m_list = (int*)calloc(config.num_m_even + config.num_m_odd, sizeof(int));
        int* m_endBatch = (int*)calloc(config.num_m_even+config.num_m_odd, sizeof(int));
	   
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
            int numRHS = temp_pointers[i] - temp_pointers[i-1];
           if (numRHS != 0)
           {
              start += numRHS;
              m_endBatch[iter] = 2 * start;
              m_list[iter] = i - 1;
              //printf("%d %d %d \n", m_list[iter], m_endBatch[iter], iter);
              iter++;
           }
        }
        config.m_list = m_list;
	    config.m_endBatch = m_endBatch;

        config.m_even_list = m_even;
        config.m_even_endBatch = m_even_endBatch;
        config.m_odd_list = m_odd;
        config.m_odd_endBatch = m_odd_endBatch;
        start = 0;
        iter = 0;
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
            int numRHS = temp_pointers[i] - temp_pointers[i-1];
            if ((((i - 1) % 2) == 0) && (numRHS!=0))
            {
                start += numRHS;
                m_even[iter] = i-1;
                m_even_endBatch[iter] = 2*start;
                //printf("%d %d %d \n", m_even[iter], m_even_endBatch[iter], iter);
                iter++;
            }
        }
    
        //current_nCols = 0;
        start = 0;
        iter = 0;
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
            int numRHS = temp_pointers[i] - temp_pointers[i-1];
            if ((((i - 1) % 2) == 1) && (numRHS!=0))
            {
                start += numRHS;
                m_odd[iter] = i-1;
                m_odd_endBatch[iter] = 2*start;
                //printf("%d %d %d \n", m_odd[iter], m_odd_endBatch[iter], iter);
                iter++;
            }

        }
        free(temp_pointers);
        config.M = (m_even[config.num_m_even - 1] > m_odd[config.num_m_odd - 1]) ? ((m_even[config.num_m_even - 1]) / 2 + 1) * 2 : ((m_odd[config.num_m_odd - 1]) / 2 + 1) * 2;// M;
        config.M = ((config.M + 31) / 32) * 32;
        config.L = (transformDirection) ? in.dims()[0] : out.dims()[0];// 3 * M / 2;
        config.Lpadded = config.L + config.M / 2 + 8; // enough for I4, standardized for matrix reuse
        if (config.Lpadded > config.Ntheta) config.Lpadded = config.Ntheta;
        if (config.projector) {
		    config.inputBufferStride = in.dims()[0];
		    config.outputBufferStride = Ntheta;
        }
        else
        {
            config.inputBufferStride = Ntheta;
		    config.outputBufferStride = out.dims()[0];
        }
	    //appContainer.input_buffer_S = (double*)in.data();
        //appContainer.buffer_S = (double*)out.data();
       initializeParallALT(&VkGPU, config, &appContainer);
        free(m_even);
        free(m_even_endBatch);
        free(m_odd);
        free(m_odd_endBatch);
        free(m_list);
        free(m_endBatch);
    };
    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    void applyImpl(Tout& out, const Tin& in){

        if (appContainer.config.Ntheta == 0)
        {
            initImpl(out, in);
        }
      Profiler::RegionFixture<5> fix("ParallaltOp::applyImpl");

      //assert(out.size() == in.size());
      //assert(out.dims()[0] == in.dims()[0]);
      assert(out.dims()[1] == in.dims()[1]);
      assert(out.dims()[2] == in.dims()[2]);
      assert(QuICC::Cuda::isDeviceMemory(out.data()));
      assert(QuICC::Cuda::isDeviceMemory(in.data()));

		/* if (temp_buffer == 0)
      {
         cudaMalloc((void**)&temp_buffer,
            appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
               sizeof(double));
      }
      if (temp_buffer2 == 0)
      {
         cudaMalloc((void**)&temp_buffer2,
            appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
               sizeof(double));
      }*/
      
      auto& memGpu = QuICC::Memory::Pensieve<Memory::Cuda::Pool>::getInstance().getMem();
      temp_buffer = reinterpret_cast<double*>(memGpu.allocate(appContainer.config.Ntheta *
                                      (appContainer.config.sizeEvenBlock +
                                         appContainer.config.sizeOddBlock) *
                                      sizeof(double)));
      temp_buffer2 = reinterpret_cast<double*>(memGpu.allocate(appContainer.config.Ntheta *
                                      (appContainer.config.sizeEvenBlock +
                                         appContainer.config.sizeOddBlock) *
                                      sizeof(double)));
      /* double* xx = (double*)calloc(2 * appContainer.config.Ntheta *
                                         (appContainer.config.sizeEvenBlock +
                                            appContainer.config.sizeOddBlock),
            sizeof(double));
        

          cudaMemcpy(xx, in.data(),
            in.dims()[0] *
               2*50 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         if (appContainer.config.projector)
         {
            for (int j = 0; j < 20; j++)
            {
               for (int i = 0; i < 5+0*in.dims()[0]; i++)
               {
                  //printf("%.17e %.17e | ", xx[2 * i+ 2*j * in.dims()[0]], xx[2 * i + 1+ 2*j * in.dims()[0]]);
               }
               //printf("\n\n");
            }
         }*/
         /* if (1)
      {

         int evenID = 0;
         int oddID = 0;
         int even_startBatch = 0;
         int odd_startBatch =
            appContainer.config
               .m_even_endBatch[appContainer.config.num_m_even - 1];
         int current_i = 0;
         
         for (int i = 0;
            i < appContainer.config.num_m_even + appContainer.config.num_m_odd;
            i++)
         {
            if (Mseq[i] % 2)
            {
               int numBatches =
                  (oddID == 0)
                     ? appContainer.config.m_odd_endBatch[0]
                     : appContainer.config.m_odd_endBatch[oddID] -
                          appContainer.config.m_odd_endBatch[oddID - 1];
               for (int l = 0; l < (numBatches / 2); l++)
               {
                  for (int j = 0; j < appContainer.config.Ntheta; j++)
                  {

                     xx2[j + odd_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l] =
                        xx[2 * j + appContainer.config.Ntheta * 2 * (current_i+l)];
                     xx2[j + odd_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l +
                         appContainer.config.Ntheta] = xx[2 * j + 1 + appContainer.config.Ntheta * 2 * (current_i+l)];
                  }
               }
               odd_startBatch += numBatches;
               current_i += numBatches / 2;
               oddID++;
            }
            else
            {
               int numBatches =
                  (evenID == 0)
                     ? appContainer.config.m_even_endBatch[0]
                     : appContainer.config.m_even_endBatch[evenID] -
                          appContainer.config.m_even_endBatch[evenID - 1];
               for (int l = 0; l < (numBatches / 2); l++)
               {
                  for (int j = 0; j < appContainer.config.Ntheta; j++)
                  {

                     xx2[j + even_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l] =
                        xx[2 * j + appContainer.config.Ntheta * 2 * (current_i+l)];
                     xx2[j + even_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l +
                         appContainer.config.Ntheta] = xx[2 * j + 1 + appContainer.config.Ntheta * 2 * (current_i+l)];
                  }
               }
               even_startBatch += numBatches;
               current_i += numBatches / 2;
               evenID++;
            }
         }

      }

        if (appContainer.config.projector)
          for (int i = 0; i < appContainer.config.Ntheta; i++)
          {
             //printf("%.3e %.3e | ", xx2[i+ Ntheta], xx2[2 * i + 1]);
         }
      cudaMemcpy(temp_buffer, xx2,
            2 * appContainer.config.Ntheta *
               (appContainer.config.sizeEvenBlock +
                  appContainer.config.sizeOddBlock) *
               sizeof(double),
            cudaMemcpyHostToDevice);*/
            
      parallALT_launchParams launchParams;
      launchParams.input_buffer_S = (double*)in.data();
      launchParams.temp_buffer_S = temp_buffer;
      launchParams.buffer_S = temp_buffer2;
      launchParams.output_buffer_S = (double*)out.data();
     // if (testAccuracy) {
				
      launchApp_parallALT(&appContainer, &launchParams);

      memGpu.deallocate(temp_buffer, appContainer.config.Ntheta *
                                     (appContainer.config.sizeEvenBlock +
                                        appContainer.config.sizeOddBlock) *
                                     sizeof(double));
       temp_buffer = 0;
    memGpu.deallocate(temp_buffer2, appContainer.config.Ntheta *
                                    (appContainer.config.sizeEvenBlock +
                                    appContainer.config.sizeOddBlock) *
                                    sizeof(double));
     temp_buffer2 = 0;
      /* double* xx = (double*)calloc(2 * appContainer.config.Ntheta *
                                         (appContainer.config.sizeEvenBlock +
                                            appContainer.config.sizeOddBlock),
            sizeof(double));
        
 if (appContainer.config.projector)
            	printf("proj\n");
            else
            	printf("int\n");
         cudaMemcpy(xx, launchParams.input_buffer_S,
            68 *
               2*8 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         //if (appContainer.config.projector)
         {
            for (int j = 0; j < 4; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 16], xx[2 * i + 1+ 2*j * 16]);
               }
               printf("\n");
            }
         }
         printf("in\n");
      cudaMemcpy(xx, launchParams.output_buffer_S,
            68 *
               2*8 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
        // if (appContainer.config.projector)
         {
            for (int j = 0; j < 4; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 68], xx[2 * i + 1+ 2*j * 68]);
               }
               printf("\n");
            }
         }
         printf("out\n");
          if (!appContainer.config.projector)
         {
            printf("I\n");
         }
         free(xx);
         cudaDeviceSynchronize();*/
      /* cudaMemcpy(xx, out.data(),
            2 * out.dims()[0] *
               528 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
      if (appContainer.config.projector)
      {
         for (int j = 0; j < 528; j++)
            {
               for (int i = 0; i < out.dims()[0]; i++)
               {
                  //printf("%.17e %.17e | ", xx[2 * i+ 2*j * out.dims()[0]], xx[2 * i + 1+ 2*j * out.dims()[0]]);
               }
               //printf("\n");
            }
      }
      /* if (1)
      {
         int evenID = 0;
         int oddID = 0;
         int even_startBatch = 0;
         int odd_startBatch =
            appContainer.config
               .m_even_endBatch[appContainer.config.num_m_even - 1];
         int current_i = 0;
        
         for (int i = 0;
            i < appContainer.config.num_m_even + appContainer.config.num_m_odd;
            i++)
         {
            if (Mseq[i] % 2)
            {
               int numBatches =
                  (oddID == 0)
                     ? appContainer.config.m_odd_endBatch[0]
                     : appContainer.config.m_odd_endBatch[oddID] -
                          appContainer.config.m_odd_endBatch[oddID - 1];
               for (int l = 0; l < (numBatches / 2); l++)
               {
                  for (int j = 0; j < appContainer.config.Ntheta; j++)
                  {

                     xx2[2 * j + appContainer.config.Ntheta * 2 * (current_i+l)] = xx[j + odd_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l];
                     xx2[2 * j + 1 + appContainer.config.Ntheta * 2 * (current_i+l)] = xx[j + odd_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l +
                         appContainer.config.Ntheta];
                     if (!appContainer.config.projector)
                     {
                        xx2[2 * j +
                            appContainer.config.Ntheta * 2 * (current_i + l)] *=
                           2 * 3.1415926535897932384626433832795029;
                        xx2[2 * j + 1 +
                            appContainer.config.Ntheta * 2 * (current_i + l)] *=
                           2 * 3.1415926535897932384626433832795029;
                     }
                  }
               }
               odd_startBatch += numBatches;
               current_i += numBatches / 2;
               oddID++;
            }
            else
            {
               int numBatches =
                  (evenID == 0)
                     ? appContainer.config.m_even_endBatch[0]
                     : appContainer.config.m_even_endBatch[evenID] -
                          appContainer.config.m_even_endBatch[evenID - 1];
               for (int l = 0; l < (numBatches / 2); l++)
               {
                  for (int j = 0; j < appContainer.config.Ntheta; j++)
                  {

                     xx2[2 * j + appContainer.config.Ntheta * 2 * (current_i+l)] = xx[j + even_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l];
                     xx2[2 * j + 1 + appContainer.config.Ntheta * 2 * (current_i+l)] =
                        xx[j + even_startBatch * appContainer.config.Ntheta + appContainer.config.Ntheta * 2 * l +
                           appContainer.config.Ntheta];
                     if (!appContainer.config.projector)
                     {
                        xx2[2 * j +
                            appContainer.config.Ntheta * 2 * (current_i + l)] *=
                           2 * 3.1415926535897932384626433832795029;
                        xx2[2 * j + 1 +
                            appContainer.config.Ntheta * 2 * (current_i + l)] *=
                           2 * 3.1415926535897932384626433832795029;
                     }
                  }
               }
               even_startBatch += numBatches;
               current_i += numBatches / 2;
               evenID++;
            }
         }
         
      }
       cudaMemcpy(out.data(), xx2,
            2 * appContainer.config.Ntheta *
               (appContainer.config.sizeEvenBlock +
                  appContainer.config.sizeOddBlock) *
               sizeof(double),
            cudaMemcpyHostToDevice);*/
      //free(xx);
      };

    void setType(int inputType, int inputDirection, int inputConvertPackedStrided)
      {
         transformType = inputType;
         transformDirection = inputDirection;
         convertPackedStrided = inputConvertPackedStrided;
      };
 private:

    /**
    * @brief parallALT configurations
    */
    mutable parallALT_configuration config = {};

    /**
    * @brief parallALT app pointers
    */
    mutable parallALT_app appContainer = {};

    mutable double* temp_buffer = 0;

    mutable double* temp_buffer2 = 0;
    /**
    * @brief parallALT app pointers
    */
    mutable PfSolve::VkGPU VkGPU = {};

    /// @brief memory resource
    /// needs shared ptr for memory pools
    /// note, this must call the dtor last
    /// otherwise we cannot dealloc data
    /// \todo consider removing shared ptr and using singleton
    std::shared_ptr<Memory::memory_resource> _mem;

    mutable int transformType = 0;
    mutable int transformDirection = 0;
    mutable int convertPackedStrided = 0;
    /// @brief Give access to base class
    //friend BinaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin, ScaleType>;

};

} // namespace Worland_parallALT
} // namespace Transform
} // namespace QuICC
